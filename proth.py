import argparse
from gmpy2 import mpz, powmod
from time import time as tm
from poprogress import simple_progress as sp
res = ""
def proth(k, n, max_a=10, vfalse=False, vprinti=False, vnores=True):
    global res
    verbose = True  # Verbose is always True
    if k > (1 << n):
        if verbose:
            res += (f"k={k} is greater than 2^n={2**n}, skipping.\n")
        return False
    if verbose and vprinti:
        res += (f"Proth test {k}*2^{n}+1")
    p_12 = mpz(k << (n - 1))
    p_1 = mpz(p_12 << 1)
    p = mpz(p_1 + 1)
    for i in range(2, max_a + 1):
        powres = powmod(i, p_12, p)
        if powres == p_1:
            if verbose:
                res += (f"Find {i}^((p-1)/2)=-1 mod p, p ({k}*2^{n}+1) is a prime with {len(str(p))} digits.\n")
            return True
        if powres != 1:
            if verbose and vfalse:
                res += (f"p ({k}*2^{n}+1) is not a prime because find {i}^((p-1)/2) != 1 and != -1 mod p\n")
            return False
    if verbose and vnores:
        res += (f"Proth test inconclusive for {k}*2^{n}+1\n")
    return False

def main():
    global res
    parser = argparse.ArgumentParser(description="Proth Prime Test")
    parser.add_argument("--min-k", type=int, default=1, help="Minimum k value (default: 1)")
    parser.add_argument("--max-k", type=int, default=10, help="Maximum k value (default: 10)")
    parser.add_argument("n", type=int, help="The n value in k*2^n+1")
    parser.add_argument("--max_a", type=int, default=10, help="The maximum value of a to test (default: 10)")
    parser.add_argument("--vfalse", action="store_true", help="Verbose output when test fails")
    parser.add_argument("--vprinti", action="store_true", help="Verbose output including initial test information")
    parser.add_argument("--vnores", action="store_true", help="Verbose output when no result is conclusive")
    
    args = parser.parse_args()
    havep = False
    t = tm()
    for k in sp(range(args.min_k // 2 * 2 + 1, args.max_k + 1, 2)):  # Ensure k is odd
        result = proth(k, args.n, args.max_a, args.vfalse, args.vprinti, args.vnores)
        if result:
            havep = True
            #print(f"{k}*2^{args.n}+1 is a prime.")
    if not havep:
        res += ("there is no prime find.\n")
    res += ("Time: %f seconds\n"%(tm()-t))
    print(res)
if __name__ == "__main__":
    main()
