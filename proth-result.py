import argparse
def main():
    parser = argparse.ArgumentParser(description="Show proth prime result (2**N*K+1)")
    parser.add_argument("-n", type=int, required=True)
    parser.add_argument("-k", type=int, required=True)
    parser.add_argument("-f", type=str, default="", required=False)
    args = parser.parse_args()
    res = 2**(args.n)*(args.k) + 1
    if args.f != "":
        print("warning: this will overwrite your file!")
        with open(args.f, "w") as file:
            file.write(str(res))
    else:
        print(res)
if __name__ =="__main__":
    main()