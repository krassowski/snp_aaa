#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function


class OutputFormatter(object):
    def __init__(self):
        self.i = 0
        self.silent = False

    def indent(self):
        self.i += 1

    def outdent(self):
        if self.i > 0:
            self.i -= 1

    def mute(self):
        self.silent = True

    def unmute(self):
        self.silent = False

    def print(self, *args, **kwargs):
        if self.silent:
            return False
        args = list(args)
        lines = str(args[0]).split('\n')
        for line in lines:
            args[0] = '\t' * self.i + str(line)
            print(*args, **kwargs)

