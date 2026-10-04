#!/usr/bin/env python3
"""Reject accidental publication before allocating build runners."""
import os
from release import tag_version


def validate(event, repository, ref_type, ref_name, release_tag='', **selections):
    if not all(isinstance(value, bool) for value in selections.values()):
        raise ValueError('Publication selections must be booleans')
    if not any(selections.values()):
        return
    if event != 'workflow_dispatch':
        raise ValueError('Publication requires a manual workflow dispatch')
    if repository.lower() != 'hallamlab/metapathways':
        raise ValueError('Publication is restricted to hallamlab/MetaPathways')
    tag = release_tag or (ref_name if ref_type == 'tag' else '')
    tag_version(tag)  # Build/artifact verification also checks its version and commit.


if __name__ == '__main__':
    selections = {}
    for name in ('PUBLISH_ANACONDA', 'PUBLISH_QUAY', 'PUBLISH_GITHUB'):
        value = os.environ.get(name, 'false').lower()
        if value not in ('true', 'false'):
            raise SystemExit(f'{name} must be true or false')
        selections[name] = value == 'true'
    validate(os.environ.get('GITHUB_EVENT_NAME', ''), os.environ.get('GITHUB_REPOSITORY', ''),
             os.environ.get('GITHUB_REF_TYPE', ''), os.environ.get('GITHUB_REF_NAME', ''),
             os.environ.get('RELEASE_TAG_INPUT', ''), **selections)
    print('Publication selections validated; unselected destinations will not publish.')
