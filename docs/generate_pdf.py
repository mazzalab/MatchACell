#!/usr/bin/env python3
"""Render a markdown doc to a clean, print-quality PDF via markdown -> HTML -> xhtml2pdf.

Usage: python docs/generate_pdf.py docs/pipeline_guide.md docs/RNASeq_pipeline_documentation.pdf

Needs: pip install --user markdown xhtml2pdf (pure Python -- no system pandoc/
LaTeX/wkhtmltopdf needed or available on this system).

xhtml2pdf/reportlab quirks worked around below (all found empirically while
building this and the sibling script in ../../circfusions/docs/generate_pdf.py):
  - Does NOT honor CSS white-space:pre/pre-wrap for <pre> layout -- embedded
    newlines and multi-space runs silently collapse. Worked around with
    explicit <br/> tags and &nbsp; padding (hardbreak_pre_blocks).
  - Does NOT support CSS page-break-inside at all (confirmed via grepping the
    installed xhtml2pdf/parser.py -- only page-break-before/after are
    parsed). Worked around by force-breaking specific blocks (tagged
    ```diagram in the source) via page-break-before, which IS supported.
  - reportlab's paragraph wrapper only breaks lines at literal space
    characters, not at '_'/'/'/'-'. Long unbroken identifiers in table cells
    are given real breakable spaces after '_'/'/' to avoid overflowing into
    the next column (add_wrap_points_to_table_code, scoped to <td> only).
"""
import re
import sys
from pathlib import Path

import markdown
from xhtml2pdf import pisa


def hardbreak_pre_blocks(html: str) -> str:
    def repl(m):
        content = m.group(2)
        content = re.sub(r" {2,}", lambda mm: " " + "&nbsp;" * (len(mm.group(0)) - 1), content)
        content = content.replace("\n", "<br/>\n")
        return m.group(1) + content + m.group(3)

    return re.sub(
        r"(<pre>(?:<code[^>]*>)?)(.*?)((?:</code>)?</pre>)",
        repl,
        html,
        flags=re.DOTALL,
    )


def add_wrap_points_to_table_code(html: str) -> str:
    # Only touch code spans long enough to plausibly overflow a narrow table
    # cell before reaching a real space -- unconditionally breaking every
    # '_'/'/' (as the sibling circfusions script does) reads fine for a genuinely
    # long, isolated token like a file path, but turns ordinary short identifiers
    # like `count_matrix` into visually broken "count_ matrix".
    LONG_ENOUGH = 18

    def fix_cell(m):
        def add_breaks(cm):
            opening, content, closing = cm.group(1), cm.group(2), cm.group(3)
            if len(content) > LONG_ENOUGH:
                content = re.sub(r"([_/])(?=\S)", r"\1 ", content)
            return opening + content + closing

        return re.sub(r"(<code[^>]*>)(.*?)(</code>)", add_breaks, m.group(0), flags=re.DOTALL)

    return re.sub(r"<td>.*?</td>", fix_cell, html, flags=re.DOTALL)


CSS = """
@font-face {
    font-family: "LibMono";
    src: url("/usr/share/fonts/liberation-mono/LiberationMono-Regular.ttf");
}
@font-face {
    font-family: "LibMono";
    font-weight: bold;
    src: url("/usr/share/fonts/liberation-mono/LiberationMono-Bold.ttf");
}
@page {
    size: A4;
    margin: 2.1cm 1.9cm 2.3cm 1.9cm;
    @frame footer {
        -pdf-frame-content: footerContent;
        bottom: 1cm; margin-left: 1.9cm; margin-right: 1.9cm; height: 1cm;
    }
}
body {
    font-family: Helvetica, sans-serif;
    font-size: 9.3pt;
    line-height: 1.45;
    color: #1a1a1a;
}
h1 {
    font-size: 20pt;
    color: #0b3d5c;
    border-bottom: 2pt solid #0b3d5c;
    padding-bottom: 6pt;
    margin-top: 0;
    margin-bottom: 14pt;
}
h2 {
    font-size: 13.5pt;
    color: #0b3d5c;
    border-bottom: 0.75pt solid #b9c9d3;
    padding-bottom: 3pt;
    margin-top: 20pt;
    margin-bottom: 8pt;
}
h3 {
    font-size: 11pt;
    color: #0b3d5c;
    margin-top: 14pt;
    margin-bottom: 6pt;
}
p { margin: 5pt 0 8pt 0; text-align: left; }
ul, ol { margin: 4pt 0 8pt 0; padding-left: 16pt; }
li { margin-bottom: 3pt; }
strong { color: #0b3d5c; }
code {
    font-family: "LibMono", Courier, monospace;
    font-size: 8.3pt;
    background-color: #eef2f5;
    color: #7a1f1f;
    padding: 0 2pt;
}
pre {
    font-family: "LibMono", Courier, monospace;
    font-size: 7.8pt;
    background-color: #f4f6f8;
    border: 0.5pt solid #cfd8de;
    padding: 6pt 8pt;
    margin: 6pt 0 10pt 0;
    line-height: 1.35;
    white-space: pre-wrap;
}
pre code { background-color: transparent; color: #1a1a1a; padding: 0; }
table {
    border-collapse: collapse;
    width: 100%;
    margin: 8pt 0 12pt 0;
    font-size: 8.2pt;
}
th, td {
    border: 0.5pt solid #b9c9d3;
    padding: 4pt 6pt;
    text-align: left;
    vertical-align: top;
}
th { background-color: #0b3d5c; color: #ffffff; }
tr:nth-child(even) td { background-color: #f4f6f8; }
hr { border: none; border-top: 0.75pt solid #b9c9d3; margin: 16pt 0; }
a { color: #0b5c94; }
#footerContent { font-size: 7.5pt; color: #6b7a83; text-align: center; }
"""

HTML_TEMPLATE = """<!DOCTYPE html>
<html>
<head><meta charset="utf-8"><title>MatchACell pipeline</title><style>{css}</style></head>
<body>
<div id="footerContent">MatchACell pipeline &mdash; page <pdf:pagenumber/> of <pdf:pagecount/></div>
{body}
</body>
</html>
"""


def main():
    src = Path(sys.argv[1])
    dst = Path(sys.argv[2])

    text = src.read_text()
    html_body = markdown.markdown(
        text,
        extensions=["extra", "tables", "fenced_code", "toc", "sane_lists"],
    )
    html_body = add_wrap_points_to_table_code(html_body)
    html_body = hardbreak_pre_blocks(html_body)
    html_body = html_body.replace(
        '<pre><code class="language-diagram">',
        '<pre style="page-break-before: always;"><code class="language-diagram">',
    )
    full_html = HTML_TEMPLATE.format(css=CSS, body=html_body)

    with dst.open("wb") as f:
        result = pisa.CreatePDF(full_html, dest=f)

    if result.err:
        print(f"ERROR: {result.err} error(s) while generating PDF", file=sys.stderr)
        sys.exit(1)
    print(f"Wrote {dst} ({dst.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
