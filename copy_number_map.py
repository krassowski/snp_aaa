#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function


class CopyNumberMap(object):

    def __init__(self, gene_start, gene_end):
        self.gene_start = gene_start
        self.gene_end = gene_end
        self.l = []
        self.whole_gene = 0

    def add(self, start, end):

        if start <= self.gene_start and end >= self.gene_end:
            self.whole_gene += 1
            return

        print('Adding inner: ', start, end)

        l = self.l

        if not l:
            l += [[start, 0], [end, 1]]
            return

        l += [[self.gene_end, 0]]
        prev = self.gene_start
        for i in range(len(l)):
            e = l[i]
            if start > prev and start < e[0]:
                l = l[:i] + [[start, e[1]]] + l[i:]
                break
            prev = e[0]
        l = l[:-1]

        l = [[self.gene_start, 0]] + l
        prev = self.gene_end
        for i in range(len(l)-1, 1, -1):
            e = l[i]
            if end < prev and end > e[0]:
                l = l[:i+1] + [[end, e[1]]] + l[i+1:]
                break
            prev = e[0]
        l = l[1:]

        l += [[self.gene_end, 0]]
        prev = self.gene_start
        for i in range(len(l)):
            e = l[i]
            if start <= prev and end >= e[0]:
                l[i][1] += 1
            prev = e[0]
        l = l[:-1]
        self.l = l

    def get_cn_list(self):
        l = self.l
        x = []
        l = [[self.gene_start, 0]] + l
        for e in l:
            if e[0] < self.gene_end and self.gene_start < e[0]:
                x += [e]
        l = l[1:]
        return x

    def show(self):
        print(self.whole_gene)

        print(self.gene_start, self.gene_end, self.gene_end - self.gene_start)
        print(self.l)

        cn_list = self.get_cn_list()

        print(cn_list)

    def get_longest(self):
        cn_list = self.get_cn_list()
        longest = 0
        count = None

        start = self.gene_start
        for e in cn_list:
            length = e[0] - start
            if length > longest:
                longest = length
                count = e[1]
            start = e[0]

        return count + self.whole_gene

    def get_max(self):

        cn_list = self.get_cn_list()
        count = max([e[1] for e in cn_list] + [0])

        return count + self.whole_gene
