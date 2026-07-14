# /// script
# requires-python = ">=3.11"
# dependencies = ["google-analytics-admin", "google-analytics-data"]
# ///
"""Pull the RBiocBook GA4 figures used in the paper.

Reports the SITE-WIDE totals (from GA4's own de-duplicated aggregate, not a
sum of the per-page column) plus the per-page breakdown, and writes a CSV
alongside the existing export in paper/data/.

CREDENTIALS: `gcloud auth application-default login` will NOT work for this.
Google blocks the analytics.readonly scope for gcloud's default OAuth client
("application denied" at the consent screen; the real error is only visible in
gcloud's stderr). Use a service account granted Viewer on the GA4 property
(GA4 Admin -> Property Access Management), then point
GOOGLE_APPLICATION_CREDENTIALS at its key -- this script picks up whatever ADC
it finds. Supplying your own OAuth client ID also works.
"""

from __future__ import annotations

import csv
from pathlib import Path

from google.analytics.admin import AnalyticsAdminServiceClient
from google.analytics.data_v1beta import BetaAnalyticsDataClient
from google.analytics.data_v1beta.types import (
    DateRange,
    Dimension,
    Metric,
    RunReportRequest,
)

START, END = "2024-05-26", "2026-06-22"
SITE_MATCH = "seandavi.github.io"
OUT = Path(
    "/Users/davsean/Documents/git/RBiocBook/paper/data/"
    f"ga4_sitewide_{START}_to_{END}.csv"
)


def find_property() -> tuple[str, str]:
    """Return (property_resource_name, display_name) for the book's GA4 property."""
    admin = AnalyticsAdminServiceClient()
    for acct in admin.list_account_summaries():
        for p in acct.property_summaries:
            if SITE_MATCH in p.display_name:
                return p.property, p.display_name
    raise SystemExit(f"No GA4 property whose name contains {SITE_MATCH!r}")


def main() -> None:
    prop, name = find_property()
    print(f"property: {prop}  ({name})")
    client = BetaAnalyticsDataClient()
    date_range = [DateRange(start_date=START, end_date=END)]
    metrics = [
        Metric(name="screenPageViews"),
        Metric(name="activeUsers"),
        Metric(name="totalUsers"),
    ]

    # Site-wide aggregate. No dimension => GA4 de-duplicates users across pages,
    # which is exactly what summing the per-page column cannot do.
    totals = client.run_report(
        RunReportRequest(property=prop, date_ranges=date_range, metrics=metrics)
    )
    print(f"\n{'=' * 46}\nSITE-WIDE {START} -> {END}\n{'=' * 46}")
    site: dict[str, str] = {}
    for header, value in zip(totals.metric_headers, totals.rows[0].metric_values):
        site[header.name] = value.value
        print(f"  {header.name:18} = {value.value}")

    # Per-page breakdown, for the record.
    pages = client.run_report(
        RunReportRequest(
            property=prop,
            date_ranges=date_range,
            dimensions=[Dimension(name="pagePath")],
            metrics=metrics + [Metric(name="userEngagementDuration")],
            limit=500,
        )
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# GA4 property", name])
        w.writerow(["# date range", f"{START} to {END}"])
        w.writerow(
            ["# SITE-WIDE (de-duplicated by GA4, not a column sum)"]
            + [f"{k}={v}" for k, v in site.items()]
        )
        cols = ["pagePath"] + [m.name for m in pages.metric_headers]
        w.writerow(cols)
        for row in pages.rows:
            w.writerow(
                [row.dimension_values[0].value]
                + [m.value for m in row.metric_values]
            )
    print(f"\nper-page rows: {len(pages.rows)}\nwrote: {OUT}")

    sum_active = sum(int(r.metric_values[1].value) for r in pages.rows)
    print(
        f"\nsanity: sum of per-page activeUsers = {sum_active} "
        f"(over-counts; site-wide activeUsers = {site.get('activeUsers')})"
    )


if __name__ == "__main__":
    main()
