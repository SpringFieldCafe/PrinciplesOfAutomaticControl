"""Export original solution JSON blocks as an editable LaTeX manuscript."""
from pathlib import Path
import os
import re

ROOT = Path(__file__).resolve().parents[1]
import json
def escape(text):
    return ''.join({'%':r'\%','&':r'\&','#':r'\#','_':r'\_','{':r'\{','}':r'\}'}.get(c,c) for c in text)
def prose(text):
    parts=re.split(r'(?<!\\)\$(.+?)(?<!\\)\$',text)
    return ''.join('$'+x+'$' if i%2 else escape(x) for i,x in enumerate(parts))
def main():
    lines=[r'% UTF-8 original solution manuscript; supplied question scans are omitted.',
           r'\documentclass[UTF8,12pt]{ctexart}',r'\usepackage{amsmath,amssymb,graphicx}',r'\begin{document}']
    # Relative graphic paths are based on this file's destination, not the caller's cwd.
    for path in sorted((ROOT/'solutions').glob('ch*.json')):
        chapter=json.loads(path.read_text(encoding='utf-8'))
        lines += ['',r'\section*{第'+str(chapter['chapter'])+'章}']
        for solution in chapter['solutions']:
            lines += ['',r'\subsection*{'+solution['id']+'}']
            for block in solution['blocks']:
                if block['type']=='p': lines += [prose(block['text']),'']
                elif block['type']=='eq': lines += [r'\[',block['latex'],r'\]','']
                elif block['type']=='figure':
                    # Use the written-assignment root as TeX working directory.
                    lines += [r'\begin{figure}[htbp]',r'\centering',r'\includegraphics[width=0.95\linewidth]{'+block['path']+'}',r'\caption{'+prose(block['caption'])+'}',r'\end{figure}','']
    lines += [r'\end{document}','']
    (ROOT/'solutions/solutions.tex').write_text('\n'.join(lines),encoding='utf-8')
if __name__=='__main__': main()
