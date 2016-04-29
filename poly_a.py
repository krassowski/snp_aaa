# -*- coding: utf-8 -*-

from __future__ import print_function


def show_pos_with_context(seq, start, end):
    return seq[:start] + '→' + seq[start:end] + '←' + seq[end:]


def has_poly_a(seq, start, end, minimal_length=12, allowed_mismatches=1, flanking=True):

    best_match = 0
    best_coords = None

    border_size = (allowed_mismatches + 2)
    border = '-' * border_size
    seq = border + seq + border
    start, end = start + border_size, end + border_size

    # -1 and + 1 allows flanking polyA to be detected too
    for starting_pos in range(start - flanking, end + flanking):
        length = 0
        mismatches = 0
        coords = []
        for d in [-1, 1]:
            pos = starting_pos
            local_mismaches = 0

            while mismatches <= allowed_mismatches:
                # print(d, pos, seq[pos], mismatches, local_mismaches)
                match = seq[pos] == 'A'
                length += match
                mismatches += not match

                pos += d
                if match:
                    local_mismaches = 0
                else:
                    local_mismaches += 1
            else:
                pos -= d * (local_mismaches + 1)
                mismatches -= (local_mismaches)
            coords.append(pos)
        # print(coords)

        if length >= best_match:
            best_match = length
            best_coords = coords

    accepted = best_match >= minimal_length
    #print(show_pos_with_context(seq, start, end), accepted)

    #print(best_coords)
    #print(seq[best_coords[0]:best_coords[1]])
    return accepted

#Persolal genome project, 100 GP - możńa sprawdzic wspołwystępowanie

if __name__ == '__main__':
    test_sequences = [
        #'     ↓  ↓' test: RASAL2, ZCRB1, RBMK2 - dlugi region, ale nie do konca sie da opisac 12-1, sprawdzic
        'AAAAAA--------',
        '---------AAAAA',
        '------AAA-AA--',
        '---AA-AAA-----',
        '------A-AAAA--',
        '-----AA-AAA---',
        'AAAAA-A-------',
        '-------A-AAAA-',
        '-----AAA------',   # False
        '---A--AA--AAA-'    # False
    ]

    for seq in test_sequences:
        print(has_poly_a(seq, 6, 9, minimal_length=5))
