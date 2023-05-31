def gm(num_seq, a, b):

    c, u_f, u_s, seq_len = input('Введите значение параметра c (c > 1/2): '), list(), list(), len(num_seq)
    divisor = max(num_seq) + 1
    
    while float(c) < 0.5:
        c = input('Некорректное значение параметра c! Попробуйте снова: ')

    c = float(c)

    q, r = c - math.log(4), c + math.sqrt(2 * c - 1)
   
    num_seq = list(map(lambda num : num / divisor, num_seq))

    for i in range(seq_len):
        if i % 2 == 0:
            u_f.append(num_seq[i])
        else: u_s.append(num_seq[i])

    upd_seq = list()
    sub_const = 1 + math.log(4.5)

    for i in range(seq_len // 2):

        f, s = u_f[i], u_s[i]

        if f == 0 or s == 0:
            continue

        V, W, Z = c * math.log(f / (1 - f)), c * (math.e ** f), f * f * s
        R = q + r * V - W

        if R + sub_const >= 4.5 * math.log(Z) or R >= math.log(Z):
            upd_seq.append(int(a + b * W))
            continue
        
        p, frac = 1 / math.sqrt(2 * c - 1), f / (1 - f)
        frac_exp = frac ** p

        if q + p * r * math.log(frac) - c * (frac ** p) + sub_const >= 4.5 * f * f * s:
            upd_seq.append(int(a + b * c * frac_exp))

    save_to_file(upd_seq, 'gm.dat')
    visualize(upd_seq, 'Гамма-распределение')