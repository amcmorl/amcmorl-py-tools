def has_vowels(word):
    has = False
    vowels = 'aeiouy'
    for vowel in vowels:
        if word.lower().__contains__(vowel):
            has = True
    return has


def unique_words(filename):

    # read in words
    f = open(filename)
    lines = f.readlines()
    line = " ".join(lines)
    words = line.split()

    # find unique words
    u = {}
    for word in words:
        # define some word rules
        if word.isalpha() and not word[0] == '\\' and not word.isupper():
            # hash the list to get unique entries only
            u[word.lower()] = 1
    uwords = u.keys()

    # more word rules - must contain a vowel
    for i, word in enumerate(uwords):
        if not has_vowels(word):
            del uwords[i]
    uwords.sort()
    return uwords


def longest(words):
    length = 0
    for word in words:
        if len(word) > length:
            maxwords = [word]
            length = len(word)
        elif len(word) == length:
            maxwords.append(word)
    return maxwords
 
