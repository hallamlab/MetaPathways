"""Normalize multi-valued EC annotations without inventing EC assignments."""
import re


def normalize_ecs(values):
    """Return distinct EC tokens in input order; preserve provisional identifiers."""
    if values is None:
        return []
    if isinstance(values, str):
        values = [values]
    result = []
    seen = set()
    for value in values:
        for token in re.split(r'[,;|]', str(value)):
            token = token.strip()
            if token and token not in seen:
                seen.add(token)
                result.append(token)
    return result
