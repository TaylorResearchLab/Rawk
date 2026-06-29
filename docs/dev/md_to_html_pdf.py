import os
import argparse
import markdown
from weasyprint import HTML



def compile_markdown(md_file_path: str,
                     pdf_file_path: str,
                     html_file_path: str):
    """
    Compiles Markdown to HTML and PDF
    """
    if not os.path.isfile(md_file_path):
        print(f"Error: The input file '{md_file_path}' does not exist.")
        return

    with open(md_file_path, "r", encoding="utf-8-sig") as f:
        md_content = f.read()

    base_dir = os.path.dirname(
        os.path.abspath(md_file_path))

    html_body = markdown.markdown(
        md_content,
        extensions=['extra', 'codehilite'],
        extension_configs={'codehilite': {'noclasses': True}}
    )

    full_html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <style>
        /* ======================================================== */
        /* GLOBAL STYLES (Applies to both HTML webpage and PDF)     */
        /* ======================================================== */
        body {{
            font-family: sans-serif;
            line-height: 1.5;
            color: #111;
            max-width: 800px; /* Optional: keeps web view easy to read */
            margin: 0 auto;   /* Centers text on wide monitors */
            padding: 20px;    /* Padding for the web view */
        }}
        pre {{
            background: #f8f8f8;
            border: 1px solid #ddd;
            padding: 12px;
            border-radius: 4px;
            white-space: pre-wrap;
            word-wrap: keep-all;
            overflow-wrap: normal;
            -webkit-line-break: after-white-space;
            line-break: strict;
            hyphens: none;
        }}
        code {{
            font-family: monospace;
            background: #f4f4f4;
        }}
        img {{
            width: auto !important;
            height: auto !important;
            max-width: 50%;
        }}

        /* ======================================================== */
        /* PRINT STYLES (WeasyPrint reads this, Web Browsers ignore)*/
        /* ======================================================== */
        @media print {{
            body {{
                max-width: 100%; /* Let text fill the full printable width */
                padding: 0;      /* Let @page margins handle spacing */
                margin: 0;
            }}
            @page {{
                size: letter;
                margin: 1in;
                @bottom-right {{
                    content: "Page " counter(page) " of " counter(pages);
                    font-family: sans-serif;
                    font-size: 9pt;
                    color: #555;
                }}
            }}
            /* Prevents a image from cutting
               in half right across a PDF page break */
            img {{
                page-break-inside: avoid;
                break-inside: avoid;
            }}
        }}
    </style>
</head>
<body>
    {html_body}

    <!-- Explicit marker to guarantee document integrity -->
    <div class="doc-end">&#9632; End of Document &#9632;</div>
</body>
</html>
"""

    # Output the static HTML file (Web browsers ignore the @media print rules)
    with open(html_file_path, "w", encoding="utf-8") as f:
        f.write(full_html)

    print(f"HTML file saved to '{html_file_path}'")

    # Output the PDF file (WeasyPrint executes the @media print rules
    # natively)
    HTML(string=full_html, base_url=base_dir).write_pdf(pdf_file_path)

    print(f"PDF file saved to '{pdf_file_path}'")



def main():
    parser = argparse.ArgumentParser(
        description="Compile Markdown to HTML and PDF."
    )
    parser.add_argument(
        "md_file",
        help="Path to the input Markdown (.md) file")
    parser.add_argument(
        "pdf_file",
        help="Path where the output PDF should be saved")
    parser.add_argument(
        "html_file",
        help="Path where the output HTML should be saved")
    
    args = parser.parse_args()
    compile_markdown(args.md_file, args.pdf_file, args.html_file)



if __name__ == "__main__":
    main()
