# /// script
# requires-python = ">=3.11"
# dependencies = ["google-analytics-admin", "google-analytics-data"]
# ///
"""Pull the RBiocBook GA4 figures used in the paper.

Reports a de-duplicated unique-reader count for THE BOOK, plus the per-page
breakdown, and writes a CSV alongside the existing export in paper/data/.

SCOPE -- read this before quoting any number this script prints.
The GA4 property "Sean Davis - web" hosts far more than the book: the BigCARE
course site, the talks site, campus-llm-kb, the agentic-AI workshop site, and
the root. As of 2026-08-01 the book was 9,098 of 17,554 property-wide views --
barely half. So the property-wide activeUsers figure is NOT a book readership
number, and an earlier version of this script reported exactly that under the
misleading label "SITE-WIDE".

The number the paper needs is unique users over pages under BOOK_PREFIX,
de-duplicated by GA4 itself. That requires BOTH:
  (a) a dimension filter on pagePath, and
  (b) NO pagePath dimension in the request -- otherwise GA4 returns one row per
      page and summing the column double-counts anyone who read two chapters.
Both queries below are labelled explicitly so the two scopes cannot be confused
again. This is the third scope error in this paper's history; hence the noise.

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
    Filter,
    FilterExpression,
    Metric,
    RunReportRequest,
)

START, END = "2024-05-26", "2026-08-01"
SITE_MATCH = "Sean Davis — web"
BOOK_PREFIX = "/RBiocBook/"
OUT = Path(
    "/Users/davsean/Documents/git/RBiocBook/paper/data/"
    f"ga4_book_{START}_to_{END}.csv"
)


def book_filter() -> FilterExpression:
    """Restrict a report to pages under BOOK_PREFIX."""
    return FilterExpression(
        filter=Filter(
            field_name="pagePath",
            string_filter=Filter.StringFilter(
                match_type=Filter.StringFilter.MatchType.BEGINS_WITH,
                value=BOOK_PREFIX,
            ),
        )
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

    # THE PAPER'S NUMBER: unique users over book pages only, de-duplicated by
    # GA4. Filtered to BOOK_PREFIX, and deliberately dimension-less so GA4
    # collapses a multi-chapter reader into one user.
    book_totals = client.run_report(
        RunReportRequest(
            property=prop,
            date_ranges=date_range,
            metrics=metrics,
            dimension_filter=book_filter(),
        )
    )
    print(f"\n{'=' * 58}\nBOOK ONLY ({BOOK_PREFIX}*)  {START} -> {END}\n{'=' * 58}")
    book: dict[str, str] = {}
    for header, value in zip(
        book_totals.metric_headers, book_totals.rows[0].metric_values
    ):
        book[header.name] = value.value
        print(f"  {header.name:18} = {value.value}")

    # CONTEXT ONLY -- the whole property, which is much more than the book.
    # Never quote this as readership.
    prop_totals = client.run_report(
        RunReportRequest(property=prop, date_ranges=date_range, metrics=metrics)
    )
    print(f"\n{'-' * 58}\nWHOLE PROPERTY (context only -- NOT book readership)\n{'-' * 58}")
    site: dict[str, str] = {}
    for header, value in zip(
        prop_totals.metric_headers, prop_totals.rows[0].metric_values
    ):
        site[header.name] = value.value
        print(f"  {header.name:18} = {value.value}")

    # Per-page breakdown for book pages, for the record.
    pages = client.run_report(
        RunReportRequest(
            property=prop,
            date_ranges=date_range,
            dimensions=[Dimension(name="pagePath")],
            metrics=metrics + [Metric(name="userEngagementDuration")],
            dimension_filter=book_filter(),
            limit=500,
        )
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# GA4 property", name])
        w.writerow(["# date range", f"{START} to {END}"])
        w.writerow(
            [f"# BOOK ONLY {BOOK_PREFIX}* (de-duplicated by GA4, not a column sum)"]
            + [f"{k}={v}" for k, v in book.items()]
        )
        w.writerow(
            ["# WHOLE PROPERTY (context only -- NOT book readership)"]
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
    landing = max(
        (int(r.metric_values[1].value) for r in pages.rows if r.dimension_values[0].value == BOOK_PREFIX),
        default=0,
    )
    print(
        "\nsanity check on the book figure -- it must sit inside these bounds:\n"
        f"  lower bound (landing page alone)   = {landing}\n"
        f"  BOOK activeUsers (GA4 de-duped)    = {book.get('activeUsers')}   <- quote this\n"
        f"  upper bound (sum of per-page)      = {sum_active}\n"
        f"  whole property, for contrast       = {site.get('activeUsers')}   <- NOT readership"
    )


if __name__ == "__main__":
    main()
