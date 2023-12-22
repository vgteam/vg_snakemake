import argparse
import fileinput

parser = argparse.ArgumentParser(description='Remove prefix from chromosome names in a SAM stream')
parser.add_argument('-f', help='fasta index listing paths and their length',
                    required=True)
parser.add_argument('-p', help='prefix to remove', required=True)
args = parser.parse_args()


def removePrefix(seqn, pref_to_rm):
    if seqn.startswith(pref_to_rm):
        return (seqn[len(pref_to_rm):])
    else:
        return (seqn)


# read paths and their length
snam = []
slen = []
path_inf = open(args.f, 'rt')
for line in path_inf:
    line = line.rstrip().split('\t')
    # # remove prefix in seq name # we shouldn't need to do this as we want the user/workfloe to prepare FASTA/FAI/DICT *without* the prefix for downstream tools to use
    # line[0] = removePrefix(line[0], args.p)
    #record length
    snam.append(line[0])
    slen.append(line[1])
path_inf.close()

# read SAM stream
write_sq = True # do we still need to write the SQ header fields
for line in fileinput.input(files='-'):
    line = line.rstrip()
    if line[0] == '@':
        # header
        if line.startswith('@SQ'):
            # we need to update those field
            if write_sq:
                # i.e. write the new seq info
                for ii in range(len(slen)):
                    print('@SQ\tSN:{}\tLN:{}'.format(snam[ii], slen[ii]))
                write_sq = False
            # and skip
            continue
        else:
            # we copy the other header fields as-is
            print(line)
    else:
        # parse line
        line = line.split('\t')
        # remove prefix
        line[2] = removePrefix(line[2], args.p)
        line[6] = removePrefix(line[6], args.p)
        # print new line
        print('\t'.join(line))
