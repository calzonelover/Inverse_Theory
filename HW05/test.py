import multiprocessing as mp

def cube(x, y):
    return x*y


if __name__ == '__main__':
    pool = mp.Pool(processes=4)
    results = [pool.apply(cube, args=(x,y)) for x in range(1,7) for y in range(7,10)]
    print(results)