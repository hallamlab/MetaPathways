"""Keep the complete CLI transcript in the selected output directory."""
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from datetime import datetime, timezone
from pathlib import Path
import shlex
import sys
import uuid


class Tee:
    def __init__(self, terminal, log):
        self.terminal, self.log = terminal, log

    def write(self, value):
        self.terminal.write(value)
        self.log.write(value)
        self.flush()
        return len(value)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

    def isatty(self):
        return False

    def fileno(self):
        return self.terminal.fileno()


@contextmanager
def transcript(directory, argv):
    folder = Path(directory).expanduser().absolute() / 'logs' / 'cli'
    folder.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    path = folder / (stamp + '-' + uuid.uuid4().hex[:8] + '.log')
    with path.open('w') as log:
        with redirect_stdout(Tee(sys.stdout, log)), redirect_stderr(Tee(sys.stderr, log)):
            print('Command: ' + shlex.join(argv))
            print(f'CLI log: {path}')
            yield
