import sys
import getopt

from calc_L_N_from_c_n import get_cn_num as cnt

def main(args):
    L, N, k, c, n = None, None, None, None, None
    try:
        opts, args = getopt.getopt(args, 'L:N:k:c:n:')
    except getopt.GetoptError as err:
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o == '-L':
            L = int(a)
        if o == '-N':
            N = int(a)
        if o == '-k':
            k = int(a)
        if o == '-c':
            c = int(a)
        if o == '-n':
            n = int(a)
    if L == None or N == None or k == None or c == None or n == None:
        print >> sys.stderr, "Missing options"
        sys.exit(2)
         
    print cnt(L, N, c, n, k)

if __name__ == '__main__':
    main(sys.argv[1:])
