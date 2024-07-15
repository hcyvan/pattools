from abc import ABC, abstractmethod
from typing import List


class Cmd(ABC):
    def __init__(self, name, help_msg):
        self.name = name
        self.help = help_msg
        self._parser = None

    def add_subparser(self, subparser):
        _parser = subparser.add_parser(self.name, help=self.help)
        self.add_argument(_parser)

    @abstractmethod
    def add_argument(self, parser):
        pass

    @abstractmethod
    def do(self, args):
        pass


class _CmdFactory:
    cmdList: List[Cmd] = []

    def __init__(self):
        pass

    def set_subparser(self, subparser):
        for _cmd in self.cmdList:
            _cmd.add_subparser(subparser)

    def run(self, args):
        for _cmd in self.cmdList:
            if _cmd.name == args.sub:
                _cmd.do(args)


CmdFactory = _CmdFactory()


def command(name, help_msg=""):
    def _cmd(cls):
        CmdFactory.cmdList.append(cls(name, help_msg))
        return cls

    return _cmd
