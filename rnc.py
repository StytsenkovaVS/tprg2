import sys
import math
import matplotlib.pyplot as plt
import numpy as np


def get_input():

    params = {}
    args = sys.argv[1:]

    for i in range(len(args)):
        flag = i == 0 or i == 1
        params[args[i][1:(2 if flag else 3)]] = args[i][(3 if flag else 4):]

    keys = params.keys()

    if 'f' not in keys:
        params['f'] = 'default.dat'

    if 'd' not in keys or 'p1' not in keys or 'p2' not in keys:
        params['d'], params['p1'], params['p2'] = 'st', 0, 100
    else: params['p1'], params['p2'] = float(params['p1']), float(params['p2'])

    return params


def save_to_file(num_seq, file_name='default.dat'):

   penult = len(num_seq) - 1

   with open(file_name, 'w') as file:
       for i in range(penult):
           file.write(str(num_seq[i]) + ',')
       file.write(str(num_seq[-1]) + '\n')


def visualize(num_seq, distrib_str, bin_num=50, cont=True):

    num_seq = np.array(num_seq)

    if cont:
        plt.hist(num_seq, bins=bin_num)
    else: plt.hist(num_seq, bins=bin_num, histtype="step")

    plt.xlabel('Значение')
    plt.ylabel('Частота')
    plt.title(distrib_str)

    plt.show()


def bi(num_seq, n, p):

    def inverse_function_method(n, p, val, bins=[], gen_bins=False):

        k, accum = 0, 0
        l_p, r_p = 1, (1 - p) ** n

        if gen_bins:

            bins = list()
            while accum < val:
                
                bin_coef = math.comb(n, k)
                accum += bin_coef * l_p * r_p
                k, l_p, r_p = k + 1, l_p * p, r_p / (1 - p)
                bins.append(bin_coef)

            return bins

        else:
            while accum < val:
                accum += bins[k] * l_p * r_p
                k, l_p, r_p = k + 1, l_p * p, r_p / (1 - p)

            return k

    divisor, n = max(num_seq) + 1, int(n)
    num_seq = list(map(lambda _num : _num / divisor, num_seq))
    num_seq = [np.random.uniform() if num_seq[i] == 0.0 else num_seq[i] for i in range(len(num_seq))]

    uni_max = max(num_seq)
    bins = inverse_function_method(n, p, uni_max, gen_bins=True)
    num_seq = list(map(lambda _num : inverse_function_method(n, p, _num, bins), num_seq))

    save_to_file(num_seq, 'bi.dat')
    visualize(num_seq, 'Биномиальное распределение', cont=False)


def ls(num_seq, a, b):

    divisor, num_seq_norm = max(num_seq) + 1, num_seq
    num_seq = list(map(lambda _num : _num / divisor, num_seq))
    num_seq = [np.random.uniform() if num_seq[i] == 0.0 else num_seq[i] for i in range(len(num_seq))]
    num_seq = list(map(lambda _num : int(a + b * math.log(_num / (1 - _num))), num_seq))
    num_seq_norm = nr(num_seq_norm, 0, 1, False)

    save_to_file(num_seq, 'ls.dat')
    visualize(num_seq, 'Логистическое распределение', 100)
    visualize(num_seq_norm, 'Нормальное распределение', 100)


def ln(num_seq, a, b):
    
    num_seq = nr(num_seq, 0, 1, False)
    num_seq = list(map(lambda _num : int(a + math.e ** (b - _num)), num_seq))

    save_to_file(num_seq, 'ln.dat')
    visualize(num_seq, 'Логнормальное распределение', 80)
    num_seq_logged = list(map(lambda _num : math.log(_num), num_seq))
    visualize(num_seq_logged, 'Логарифм логнормального распределения', 80)


def gm(num_seq, a, b):

    c, divisor = 1, max(num_seq) + 1

    c = int(c)

    num_seq = list(map(lambda num : num / divisor, num_seq))
    num_seq = [num_seq[i:i+c] for i in range(0, len(num_seq), c)]

    if len(num_seq[-1]) != c:
        del num_seq[-1]

    for i in range(len(num_seq)):
        if 0.0 in num_seq[i]:
            num_seq[i] = [np.random.uniform() if num_seq[i][j] == 0.0 else num_seq[i][j] for j in range(c)]

    num_seq = list(map(lambda _list : int(a - b * math.log(np.prod(list(map(lambda _num : 1 - _num, _list))))), num_seq))

    save_to_file(num_seq, 'gm.dat')
    visualize(num_seq, 'Гамма-распределение', 100)


def nr(num_seq, mu, sigma, targeted=True):

    def box_muller(f, s, func):
        temp = mu + sigma * math.sqrt(-2 * math.log(1 - f)) * func(2 * math.pi * s)
        return (int(temp) if targeted else temp)

    def aux(u_f, u_s, func):
        return list(map(lambda f, s : box_muller(f, s, func), u_f, u_s))

    divisor = max(num_seq) + 1
    num_seq, u_f, u_s = list(map(lambda num : num / divisor, num_seq)), list(), list()

    for i in range(len(num_seq)):
        if i % 2 == 0:
            u_f.append(num_seq[i])
        else: u_s.append(num_seq[i])

    u_f, u_s = aux(u_f, u_s, math.cos), aux(u_f, u_s, math.sin)

    num_seq = u_f + u_s

    if targeted:
        save_to_file(num_seq, 'nr.dat')
        visualize(num_seq, 'Нормальное распределение', 100)

    return num_seq


def ex(num_seq, a, b):

    divisor = max(num_seq) + 1
    num_seq = [num_seq[i] for i in range(len(num_seq)) if num_seq[i] != 0]
    num_seq = list(map(lambda num : int(-b * math.log(num) + a), list(map(lambda num : num / divisor, num_seq))))

    save_to_file(num_seq, 'ex.dat')
    visualize(num_seq, 'Общее экспоненциальное распределение', 100)


def tr(num_seq, a, b):

   divisor = max(num_seq) + 1
   num_seq, u_f, u_s = list(map(lambda num : num / divisor, num_seq)), list(), list()

   for i in range(len(num_seq)):
       if i % 2 == 0:
           u_f.append(num_seq[i])
       else: u_s.append(num_seq[i])
   
   num_seq = list(map(lambda f, s : int(a + b * (f + s - 1)), u_f, u_s))

   save_to_file(num_seq, 'tr.dat')
   visualize(num_seq, 'Треугольное распределение', 100)


def st(num_seq, a, b):

    divisor = max(num_seq) + 1
    num_seq = list(map(lambda num : int(num / divisor * b + a), num_seq))
    
    save_to_file(num_seq, 'st.dat')
    visualize(num_seq, 'Стандартное равномерное распределение с заданным интервалом')


def main():
    
    params = get_input()
    file, distrib, p1, p2 = params.values()
   
    f = open(file, "r")
    num_seq = list(map(int, (f.read()).split(',')))

    if distrib=='st':
        st(num_seq, p1, p2)
    elif distrib=='tr':
        tr(num_seq, p1, p2)
    elif distrib=='ex':
        ex(num_seq, p1, p2)
    elif distrib=='nr':
        nr(num_seq, p1, p2)
    elif distrib=='gm':
        gm(num_seq, p1, p2)
    elif distrib=='ln':
        ln(num_seq, p1, p2)
    elif distrib=='ls':
        ls(num_seq, p1, p2)
    elif distrib=='bi':
        bi(num_seq, p1, p2)
    else:
        print('Такого распределения нет!')
        return


if __name__ == '__main__':
    main()
