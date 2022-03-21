# Strings for SI LaTeX files

header = """\\documentclass[10pt,oneside]{article}
\\usepackage[utf8]{inputenc}
\\usepackage[margin=1in]{geometry}
\\usepackage{multicol}
\\usepackage{booktabs}
\\usepackage[version=4]{mhchem}
\\usepackage{pdflscape}

\\renewcommand{\\arraystretch}{1.2}
\\extrafloats{100}

\\renewcommand{\\thepage}{S\\arabic{page}}
\\renewcommand{\\thesection}{S\\arabic{section} }
\\renewcommand{\\thetable}{S\\arabic{table}}
\\renewcommand{\\thefigure}{S\\arabic{figure}}

\\title{Supplementary Information}
\\author{}
\\date{}

\\begin{document}

\\maketitle

\\tableofcontents

\\newpage

"""

footer = """
\\end{document}
"""
