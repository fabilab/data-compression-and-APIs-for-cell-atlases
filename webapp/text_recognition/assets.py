# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/05/22
content:    Some static assets for text interpretation
'''
def get_numbers():
    numbers = {
        'one': 1,
        'two': 2,
        'three': 3,
        'four': 4,
        'five': 5,
        'six': 6,
        'seven': 7,
        'eight': 8,
        'nine': 9,
        'ten': 10,
        'eleven': 11,
        'twelve': 12,
        'thirteen': 13,
        'fourteen': 14,
        'fifteen': 15,
        'sixteen': 16,
        'seventeen': 17,
        'eighteen': 18,
        'nineteen': 19,
        'twenty': 20,
        'thirty': 30,
        'fourty': 40,
        'fifty': 50,
        'sixty': 60,
        'seventy': 70,
        'eighty': 80,
        'ninety': 90,
        'hundred': 100,
    }
    numbers_inv = {val: key for key, val in numbers.items()}
    for n in range(21, 100):
        if n in numbers_inv:
            continue
        decade = n - n % 10
        digit = n % 10
        ntext = numbers_inv[decade] + numbers_inv[digit]
        numbers[ntext] = n

    numbers = list(numbers.items())
    numbers.sort(key=lambda x: x[1])

    return numbers


numbers = get_numbers()
