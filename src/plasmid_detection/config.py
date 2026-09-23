"""Load config/config.yaml, resolve ${dotted.key} references, and expose it to shell stages.

python -m plasmid_detection.config shell         # export KEY='value' lines for eval
python -m plasmid_detection.config get db.gtdb.dir
"""

import argparse
import re
import shlex
import sys
from pathlib import Path

import yaml

DEFAULT_CONFIG = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
_REF = re.compile(r"\$\{([A-Za-z0-9_.]+)\}")


def _lookup(tree, dotted):
    node = tree
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            raise KeyError(f"unknown config key: {dotted}")
        node = node[part]
    return node


def _resolve(tree, value, seen=()):
    if isinstance(value, dict):
        return {k: _resolve(tree, v, seen) for k, v in value.items()}
    if not isinstance(value, str):
        return value

    def sub(match):
        key = match.group(1)
        if key in seen:
            raise ValueError(f"circular config reference: {' -> '.join(seen + (key,))}")
        target = _lookup(tree, key)
        if isinstance(target, dict):
            raise ValueError(f"config reference to a section, not a value: {key}")
        return str(_resolve(tree, target, seen + (key,)))

    return _REF.sub(sub, value)


def load(path=DEFAULT_CONFIG):
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    return _resolve(raw, raw)


def flatten(tree, prefix=""):
    for key, value in tree.items():
        name = f"{prefix}_{key}" if prefix else key
        if isinstance(value, dict):
            yield from flatten(value, name)
        else:
            yield name.upper(), value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("shell", help="print export statements for every key")
    get = sub.add_parser("get", help="print one resolved value")
    get.add_argument("key")
    args = parser.parse_args(argv)

    cfg = load(args.config)
    if args.cmd == "shell":
        for name, value in flatten(cfg):
            print(f"export {name}={shlex.quote(str(value))}")
    else:
        print(_lookup(cfg, args.key))
    return 0


if __name__ == "__main__":
    sys.exit(main())
