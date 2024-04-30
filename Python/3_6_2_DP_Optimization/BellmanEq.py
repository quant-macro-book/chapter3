def BellmanEq(m,capital,kprime,v_interp):
    """
    -----------------------------------------------
    === (k,k')を所与とした時にベルマン方程式を返す関数 ===
    -----------------------------------------------
    <input>
    ・m: パラメータ等を格納したコンストラクタ
    ・capital: 今期の資本保有量(k)
    ・kprime: 来期の資本保有量(k')
    ・v_: 補間した来期の価値関数
    (※Matlabコード上ではBellmanEq関数内で補間した価値関数を定義しているが、
    今回は明示的に関数の引数として扱う(capitalも同様))
    <output>
    ・value: (k,k')に対するベルマン方程式
    (※コード上では"最小化問題"を解くので符号を反転させた値を返す)
    """
    from CRRA import CRRA

    alpha, beta = m.alpha, m.beta

    if kprime < 0: #(1): k'は正の値しか取らないので、ペナルティを与えてその値が選ばれないようにする
        
        value = -1000000.0
    
    else:

        wealth = capital ** alpha #今期の資本保有量で実現される生産量
        cons = wealth - kprime #消費量
        util = CRRA(m,cons) #効用水準
        value = util + beta*v_interp(kprime) #u(c)+βV_{t+1}(k')

    value = -1*value #(2): コード上では"最小化問題"を解くので符号を反転させた値を返す

    return value

    

    



    