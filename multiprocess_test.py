# -*- coding: utf-8 -*-


import time
from multiprocessing import Pool

#並列処理させる関数
def nijou(x):
    print('input: %d' % x)
    time.sleep(2)
    retValue = x * x
    print('double: %d' % (retValue))
    return(retValue)

if __name__ == "__main__":
    p = Pool(4) # プロセス数を4に設定
    result = p.map(nijou, range(10))  # nijou()に0,1,..,9を与えて並列演算
    print(result)