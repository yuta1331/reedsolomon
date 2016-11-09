#!/usr/bin/env python

def _return_width_list(num):
	i = 1
	while True:
		if ((i-1)**2 < num) & (num <= i**2):
			if num % i:
				return [i, (num // i) + 1]
			else:
				return [i, num // i]
		if i == 1000:
			print('ERROR:subset._return_width: number of bits is too big.')
			return 'ERROR:subset._return_width: number of bits is too big.'
		i += 1

def interleaver(lists, option):
	width = _return_width_list(len(lists))

	if option == 'interleave':
		while len(lists) % width[0]:
			lists.append(None)
	elif option == 'deinterleave':
		if len(lists) % width[1]:
			tmp = list()
			cnt_enough = len(lists) % width[0]
			cnt_not_enough = 0
			for cnt, l in enumerate(lists):
				tmp.append(l)
				if cnt_enough:
					if not (cnt+1) % width[1]:
						cnt_enough -= 1
				else:
					cnt_not_enough += 1
					if cnt_not_enough == width[1] - 1:
						tmp.append(None)
						cnt_not_enough = 0
			lists = tmp
	else:
		print('ERROR:subset.interleaver: option is wrong.')
		return 'ERROR:subset.interleaver: option is wrong.'

	result = list()
	width = width[option=='deinterleave']
	for i in range(width):
		for j in range(len(lists) // width):
			if not lists[width*j + i] == None:
				result.append(lists[width*j + i])
	return result
########## Usage ##########
#lists = list(range(18))  # [0,1,...,17]
#lists = interleaver(lists, 'interleave')  # [0,5,...,9,14]
#lists = interleaver(lists, 'deinterleave')  # [0,1,...,17]

########## reedsolomon(refering to reedsolo.py: library) ##########
#str_input is bynary data
#reedsolomon(str_input)

class ReedSolomonError(Exception):
    pass

gf_exp = [1] * 512
gf_log = [0] * 256
x = 1
for i in range(1, 255):
    x <<= 1
    if x & 0x100:
        x ^= 0x11d
    gf_exp[i] = x
    gf_log[x] = i
for i in range(255, 512):
    gf_exp[i] = gf_exp[i - 255]

def gf_mul(x, y): # x * y on gf_exp
    if x == 0 or y == 0:
        return 0
    return gf_exp[gf_log[x] + gf_log[y]]

def gf_div(x, y): # x / y on gf_exp
    if y == 0:
        raise ZeroDivisionError()
    if x == 0:
        return 0
    return gf_exp[gf_log[x] + 255 - gf_log[y]]

def gf_poly_scale(p, x): # p_list * x on gf_exp # [(p[0] * x), (p[1] * x), ...]
    return [gf_mul(p[i], x) for i in range(0, len(p))]

def gf_poly_add(p, q):
    r = [0] * max(len(p), len(q))
    for i in range(0, len(p)):
        r[i + len(r) - len(p)] = p[i]  # or r[i]...?
    for i in range(0, len(q)):
        r[i + len(r) - len(q)] ^= q[i] # or r[i]...?
    return r
# p       1011, 0011        p 1011, 0011
# q 1001, 1010, 0101   or   q 1001, 1010, 0101
# r 1001, 0001, 0110        r 0010, 1001, 0101

def gf_poly_mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for j in range(0, len(q)):
        for i in range(0, len(p)):
            r[i + j] ^= gf_mul(p[i], q[j])
    return r

def gf_poly_eval(p, x):
    y = p[0]
    for i in range(1, len(p)):
        y = gf_mul(y, x) ^ p[i]
    return y

########## reedsolomon ##########
# nsym is the number of bytes in the error correction code.
# the algorithm can correct up to ``nsym/2`` of the errors in the message.

def rs_generator_poly(nsym):
    g = [1]
    for i in range(0, nsym):
        g = gf_poly_mul(g, [1, gf_exp[i]])
    return g

def rs_encode_msg(msg_in, nsym): ####### for encoding message #######
    if len(msg_in) + nsym > 255:
        raise ValueError("message too long")
    gen = rs_generator_poly(nsym)
    msg_out = [0] * (len(msg_in) + nsym)
    msg_out[:len(msg_in)] = msg_in
    for i in range(0, len(msg_in)):
        coef = msg_out[i]
        if coef != 0:
            for j in range(0, len(gen)):
                msg_out[i + j] ^= gf_mul(gen[j], coef)
    msg_out[:len(msg_in)] = msg_in
    return msg_out

def rs_calc_syndromes(msg, nsym):
    return [gf_poly_eval(msg, gf_exp[i]) for i in range(nsym)]

def rs_correct_errata(msg, synd, pos):
    # calculate error locator polynomial
    q = [1]
    for i in range(0, len(pos)):
        x = gf_exp[len(msg) - 1 - pos[i]]
        q = gf_poly_mul(q, [x, 1])
    # calculate error evaluator polynomial
    p = synd[0:len(pos)]
    p.reverse()
    p = gf_poly_mul(p, q)
    p = p[len(p) - len(pos):len(p)]
    # formal derivative of error locator eliminates even terms
    q = q[len(q) & 1:len(q):2]
    # compute corrections
    for i in range(0, len(pos)):
        x = gf_exp[pos[i] + 256 - len(msg)]
        y = gf_poly_eval(p, x)
        z = gf_poly_eval(q, gf_mul(x, x))
        msg[pos[i]] ^= gf_div(y, gf_mul(x, z))

def rs_find_errors(synd, nmess):
    # find error locator polynomial with Berlekamp-Massey algorithm
    err_poly = [1]
    old_poly = [1]
    for i in range(0, len(synd)):
        old_poly.append(0)
        delta = synd[i]
        for j in range(1, len(err_poly)):
            delta ^= gf_mul(err_poly[len(err_poly) - 1 - j], synd[i - j])
        if delta != 0:
            if len(old_poly) > len(err_poly):
                new_poly = gf_poly_scale(old_poly, delta)
                old_poly = gf_poly_scale(err_poly, gf_div(1, delta))
                err_poly = new_poly
            err_poly = gf_poly_add(err_poly, gf_poly_scale(old_poly, delta))
    errs = len(err_poly) - 1
    if errs * 2 > len(synd):
        raise ReedSolomonError("Too many errors to correct")
    # find zeros of error polynomial
    err_pos = []
    for i in range(0, nmess):
        if gf_poly_eval(err_poly, gf_exp[255 - i]) == 0:
            err_pos.append(nmess - 1 - i)
    if len(err_pos) != errs:
        return None    # couldn't find error locations
    return err_pos

def rs_forney_syndromes(synd, pos, nmess):
    fsynd = list(synd)      # make a copy
    for i in range(0, len(pos)):
        x = gf_exp[nmess - 1 - pos[i]]
        for i in range(0, len(fsynd) - 1):
            fsynd[i] = gf_mul(fsynd[i], x) ^ fsynd[i + 1]
        fsynd.pop()
    return fsynd

def rs_correct_msg(msg_in, nsym): ####### for correcting message #######
    if len(msg_in) > 255:
        raise ValueError("message too long")
    msg_out = list(msg_in)     # copy of message
    # find erasures
    erase_pos = []
    for i in range(0, len(msg_out)):
        if msg_out[i] < 0:
            msg_out[i] = 0
            erase_pos.append(i)
    if len(erase_pos) > nsym:
        raise ReedSolomonError("Too many erasures to correct")
    synd = rs_calc_syndromes(msg_out, nsym)
    if max(synd) == 0:
        return msg_out[:-nsym]  # no errors
    fsynd = rs_forney_syndromes(synd, erase_pos, len(msg_out))
    err_pos = rs_find_errors(fsynd, len(msg_out))
    if err_pos is None:
        raise ReedSolomonError("Could not locate error")
    rs_correct_errata(msg_out, synd, erase_pos + err_pos)
    synd = rs_calc_syndromes(msg_out, nsym)
    if max(synd) > 0:
        raise ReedSolomonError("Could not correct message")
    return msg_out[:-nsym]

# binary string -> decimal list
def BtoD(str_in):
    deci = []
    for i in range(0, len(str_in), 8):
        deci.append(int(str_in[i:i+8], 2))
    return deci

# decimal list -> binary string
def DtoB(list_in):
    bina = str()
    for i in range(len(list_in)):
        list_in[i] = format(list_in[i], 'b')
        while len(list_in[i]) < 8:
            list_in[i] = '0' + list_in[i]
        bina += list_in[i]
    return bina

########## API ##########

class RSCodec(object):
    
    def __init__(self, nsym=10):
        self.nsym = nsym

    def encode(self, data):
        return DtoB(rs_encode_msg(BtoD(data), self.nsym))
        
    def decode(self, data):
        return DtoB(rs_correct_msg(BtoD(data), self.nsym))
        
