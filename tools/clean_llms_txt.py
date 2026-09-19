"""Post-render tidy-up of the llms.txt that Quarto generates for the book.

Quarto copies book chapter titles into llms.txt with their Pandoc span markup
(``[[3]{.chapter-number}  [R mechanics]{.chapter-title}]``). This rewrites them
as plain titles ("3 R mechanics") and adds the one-line summary blockquote the
llms.txt format expects under the H1. It is a no-op when llms.txt is absent
(e.g. a pdf- or epub-only render).
"""

import os
import re
from pathlib import Path

SUMMARY = (
    "> An open textbook introducing R, statistics, and Bioconductor for "
    "biologists. Each link is a Markdown version of one chapter, with its code."
)

SPAN = re.compile(r"\[([^\[\]]*)\]\{\.chapter-(?:number|title)\}")


def clean(text: str) -> str:
    lines = []
    for line in text.splitlines():
        # [[3]{.chapter-number}  [R mechanics]{.chapter-title}](url) -> [3 R mechanics](url)
        line = SPAN.sub(r"\1", line)
        line = re.sub(r"(?<=\[)(\d+|[A-Z])\s{2,}", r"\1 ", line)
        lines.append(line)
    if SUMMARY not in lines:
        h1 = next((i for i, l in enumerate(lines) if l.startswith("# ")), None)
        if h1 is not None:
            lines[h1 + 1:h1 + 1] = ["", SUMMARY]
    return "\n".join(lines) + "\n"


def main() -> None:
    out_dir = Path(os.environ.get("QUARTO_PROJECT_OUTPUT_DIR", "_book"))
    llms = out_dir / "llms.txt"
    if llms.exists():
        llms.write_text(clean(llms.read_text(encoding="utf-8")), encoding="utf-8")


if __name__ == "__main__":
    main()
