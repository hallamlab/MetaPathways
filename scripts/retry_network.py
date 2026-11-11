#!/usr/bin/env python3
"""Retry transient network failures without overwriting remote artifacts."""
import re
import subprocess
import sys
import tempfile
import time


def transient(output):
    # Authentication and existing-artifact conflicts require operator review.
    if re.search(r'\b(?:401|403|409)\b|already exists|Unauthorized|Forbidden', output, re.I):
        return False
    return bool(re.search(
        r'\b(?:429|500|502|503|504|520|521|522|523|524)\b|'
        r'Temporary failure resolving|Temporary failure in name resolution|'
        r'Unable to find the server|Could not resolve|Connection reset|'
        r'Connection timed out|Read timed out', output, re.I))


def run(command):
    for attempt in range(3):
        # Keep subprocess output visible while retaining it for classification.
        with tempfile.TemporaryFile(mode='w+') as log:
            with subprocess.Popen(command, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, text=True) as process:
                for line in process.stdout:
                    print(line, end='', flush=True)
                    log.write(line)
            log.seek(0)
            output = log.read()
        if process.returncode == 0:
            return
        if attempt == 2 or not transient(output):
            raise subprocess.CalledProcessError(process.returncode, command)
        delay = 15 * (attempt + 1)
        print(f'Transient network failure; retry {attempt + 2}/3 in {delay}s.', flush=True)
        time.sleep(delay)


if __name__ == '__main__':
    try:
        run(sys.argv[1:])
    except subprocess.CalledProcessError as error:
        sys.exit(error.returncode)
