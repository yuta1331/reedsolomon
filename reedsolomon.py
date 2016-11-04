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


