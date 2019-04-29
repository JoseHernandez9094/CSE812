import pandas as pd


def main():
    des = "Data/DIM3__SEED__"
    fast,name = [],[]
    solve = 0

    for i in range(1, 101):
        f = des + str(i) + '/data.csv'
        low = 201

        df = pd.read_csv(f, sep=',')
        for r in df.itertuples():
            cur = int(r._2)
            if cur < low:
                low = cur

        if low < 201:
            solve += 1
            fast.append(low)
            name.append('3x3 Grid')

    print('Cnt: ', solve)

    df = pd.DataFrame({'i':fast, 'Instruction':name})
    df.to_csv('fast.csv', index=False)

if __name__ == "__main__":
    main()