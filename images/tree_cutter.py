"""
takes latex document and associated file tree, replaces all \\input{file} with the contents of the tex file
also optionally deletes all comments
"""

import argparse, os

PARSER = argparse.ArgumentParser()
PARSER.add_argument('--main', action='store', required=True, help='main tex file')
PARSER.add_argument('--dir', action='store', required=False,
                    help='root of file tree to search (defaults to parent dir of main tex file)')
PARSER.add_argument('--output-file', action='store', required=False,
                    help='where to save new tex file (defaults to /output/<main tex file name>)')
PARSER.add_argument('--comments-remove', action='store_true', required=False,
                    help="delete anything line beginning with '%%'")
PARSER.add_argument('--cmds-remove', action='store', nargs='+', required=False, default=[],
                    help="delete anything with the corresponding command (e.g. '--rm-cmd todo' will remove any \\todo{...})")

args = PARSER.parse_args()
main_file = args.main
root_dir = args.dir if args.dir else os.path.dirname(main_file)
output_file = args.output_file if args.output_file else os.path.join('output',
                                                                     os.path.basename(main_file)
                                                                     )
if not os.path.exists(os.path.dirname(output_file)):
    os.makedirs(os.path.dirname(output_file))

f = open(main_file, 'r')
txt = f.read()
f.close()
include_cmds = ['\\input{', '\\include{']
things = include_cmds + ['\\' + cmd + '{' for cmd in args.cmds_remove]
for thing in things:
    while thing in txt:
        bgn = txt.index(thing)
        q = txt[bgn + len(thing):]
        j = -1
        c = 1
        while c > 0:
            j += 1
            if not q[:j].endswith('\\'):
                if q[j] == '}':
                    c -= 1
                if q[j] == '{':
                    c += 1
        q = q[:j]
        if thing in include_cmds:
            fn = os.path.join(root_dir, *os.path.split(q)) + '.tex'
            print('extracting:', fn)
            if not os.path.exists(fn):
                print('NOT FOUND:', fn)
                txt = txt[:bgn] + txt[bgn + len(thing) + len(q) + 1:]
                continue
            f = open(fn, 'r')
            txt = txt[:bgn] + f.read() + txt[bgn + len(thing) + len(q) + 1:]
            f.close()
        else:
            txt = txt[:bgn] + txt[bgn + len(thing) + len(q) + 1:]

if args.comments_remove:
    i = 0
    while i < len(txt):
        if txt[i] == '%' and not txt[:i].endswith('\\'):
            lgth = txt[i:].index('\n') + 1
            txt = txt[:i] + txt[i + lgth:]
        else:
            i += 1

f = open(output_file, 'w')
f.write(txt)
f.close()
print('saved to', output_file)
