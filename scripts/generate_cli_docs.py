#!/usr/bin/env python3
"""Generate the GitHub CLI reference from the installed source parser help."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
COMMANDS=('prepare_test','run','analysis_wf','build_db','mag_split','build_pt','ptools','report')


def stable_action_invocation(self, action):
    """Keep reference help in Python 3.11 form across argparse versions.

    Python 3.13 elides repeated alias metavars. Limit this compatibility shim
    to the documentation subprocess; the installed CLI keeps native formatting.
    """
    if not action.option_strings:
        default = self._get_default_metavar_for_positional(action)
        return self._metavar_formatter(action, default)(1)[0]
    if action.nargs == 0:
        return ', '.join(action.option_strings)
    default = self._get_default_metavar_for_optional(action)
    args = self._format_args(action, default)
    return ', '.join(option + ' ' + args for option in action.option_strings)


def render(maintainer=False):
    env=dict(os.environ, PYTHONPATH=str(ROOT), COLUMNS='100')
    title = 'Maintainer CLI reference' if maintainer else 'Complete user CLI reference'
    guide = '[Development](development.md) · [Reaction screening](ptools-maintenance.md)' if maintainer else '[User guide](index.md) · [Advanced user guide](advanced.md)'
    text=f'# {title}\n\n{guide}\n\n'
    text+='Run `metapathways COMMAND --help` for the help of your installed revision.\n\n'
    for command in (('screen_pt',) if maintainer else COMMANDS):
        code=('import sys, argparse; '
              'from scripts.generate_cli_docs import stable_action_invocation; '
              'argparse.HelpFormatter._format_action_invocation=stable_action_invocation; '
              'from metapathways.pipeline import main; sys.argv=["metapathways"]+sys.argv[1:]; main()')
        result=subprocess.run([sys.executable,'-c',code,command,'--help'],cwd=ROOT,env=env,check=True,text=True,capture_output=True)
        text+=f'## {command}\n\n```text\n{result.stdout.rstrip()}\n```\n\n'
    return '\n'.join(line.rstrip() for line in text.splitlines()).rstrip() + '\n'


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    for maintainer, name in ((False, 'cli-reference.md'), (True, 'developer-cli-reference.md')):
        target=ROOT/'docs'/name; expected=render(maintainer)
        if args.check:
            if not target.is_file() or target.read_text()!=expected:
                raise SystemExit(f'{name} is stale; run python scripts/generate_cli_docs.py')
            print(f'{name} matches current parser help.')
        else:
            target.write_text(expected)

if __name__=='__main__':main()
