'''
Amr Yehia Zakaria 20198062
Mohamed Ibrahim Shawky 20198068
Mahmoud Khaled Helmy 20188045
Yomna Farid Mohamed 20198102 
Lamiaa Karam Zaki 20198120
Ghada Mamdouh Ahmed 20198063
'''

import random
import numpy as np
random.seed(1)

x = int(input('Median string press 1, PSSM press 2 : '))
if x == 1:
    def get_sequences():
        num = int(input("Press 1 to use the rawDNA file or 2 to generate random sequences:"))
        if num == 1:
            file_path = "rawDNA.txt"
            sequences, t, n, l = read_file_med(file_path)
            return sequences, t, n, l
        elif num == 2:
            sequences, t, n, l = generate_random_seq_med()
            return sequences, t, n, l
        else:
            print("Invalid input. Please enter either 1 or 2.")
            return get_sequences()


    def read_file_med(path):
        try:
            f = open(path, "r")
            data = []
            for i in f:
                i = i.split()
                data.append(i)

            # close the file
            f.close()
            dnas = data[1:]
            dnas = [gene[0] for gene in dnas]
            dnas = [n.upper() for n in dnas]
            number_of_dnas = int(data[0][0])
            len_of_dnas = int(data[0][1])
            len_of_patterns = int(data[0][2])
            return dnas, number_of_dnas, len_of_dnas, len_of_patterns
        except FileNotFoundError:
            print('File not found, Get the random multiple sequences')
            sequences, t, n, l = generate_random_seq()
            return sequences, t, n, l


    def generate_random_seq_med():
        len_of_patterns = int(input("Enter the length of Patterns: "))
        number_of_dnas = int(input("Enter number of sequences:"))
        len_of_dnas = int(input("Enter number of nucleotides:"))
        nucleotides = ['A', 'C', 'G', 'T']
        sequences = [''.join(random.choice(nucleotides) for _ in range(len_of_dnas)) for _ in range(number_of_dnas)]
        return sequences, number_of_dnas, len_of_dnas, len_of_patterns


    CONVERT_STR = ['A', 'C', 'G', 'T']


    def create_empy_string(word):
        result = [' '] * len(word)
        for i in range(len(word)):
            result[i] = CONVERT_STR[word[i]]
        return ''.join(result)


    def bypass(short_seq, level, max_num):
        i = len(short_seq) - 1
        while i >= 0:
            if short_seq[i] < max_num:
                short_seq[i] += 1
                level = i + 1
                return level
            i -= 1
        level = 0
        return level


    def next_vertex(short_seq, level, max_num):
        if level < len(short_seq):
            short_seq[level] = 0
            level += 1
            return level
        else:
            level = bypass(short_seq, level, max_num)
            return level


    def Distance(seq, word, level):
        min_distance = float('inf')
        for i in range(len(seq) - len(word) + 1):
            distance = 0
            for j in range(level):
                if seq[i + j] != word[j]:
                    distance += 1
            if distance < min_distance:
                min_distance = distance
        return min_distance


    def TotalDistance(dnas, word, level):
        total_distance = 0
        for i in range(len(dnas)):
            print('Sequence: {}, index: {}, Word: {}'.format(dnas[i], i, word))
            total_distance += Distance(dnas[i], word, level)
        return total_distance


    def branch_and_bound_median_string_search(dnas, len_of_patterns):
        word = [0] * len_of_patterns
        best_word = [0] * len_of_patterns
        best_distance = float('inf')
        level = 1
        optimistic_distance = 0
        max_word_number = 3
        total_distance = 0
        while level > 0:
            if level < len_of_patterns:
                prefix = create_empy_string(word)
                optimistic_distance = TotalDistance(dnas, prefix, level)
                if optimistic_distance > best_distance:
                    level = bypass(word, level, max_word_number)
                else:
                    level = next_vertex(word, level, max_word_number)
            else:
                prefix = create_empy_string(word)
                total_distance = TotalDistance(dnas, prefix, level)
                if total_distance < best_distance:
                    best_distance = total_distance
                    best_word = word.copy()
                level = next_vertex(word, level, max_word_number)
        return create_empy_string(best_word)

    def main():
        dnas, number_of_dnas, len_of_dnas, len_of_patterns = get_sequences()
        best_word = branch_and_bound_median_string_search(dnas, len_of_patterns)
        print("The string:", best_word)

    main()
elif x == 2:
    def get_multi_aligned_seq_mat():
        num = int(input("Press 1 to use the PSSMData file or 2 to generate random sequences:"))
        if num == 1:
            file_path = "PSSMData.txt"
            multi_aligned_seq_matrix, t, n = read_file(file_path)
            return multi_aligned_seq_matrix, t, n
        elif num == 2:
            multi_aligned_seq_matrix, t, n = generate_random_seq()
            return multi_aligned_seq_matrix, t, n
        else:
            print("Invalid input. Please enter either 1 or 2.")
            return get_multi_aligned_seq_mat()


    '''read the PSSMData.txt file'''


    def read_file(path):
        try:
            file = open(path, "r")
            words = []
            for line in file:
                # split the line into words
                temp = line.split()
                words.append(temp)

            # close the file
            file.close()
            # Extract the sequences from the list
            sequences = words[1:]
            # Convert the sequences to a matrix
            multi_aligned_seq = np.array([list(seq[0]) for seq in sequences])
            return multi_aligned_seq, int(words[0][0]), int(words[0][1])
        except FileNotFoundError:
            print('File not found, Get the random multiple sequences')
            multi_aligned_seq_matrix, t, n = generate_random_seq()
            return multi_aligned_seq_matrix, t, n


    def generate_random_seq():
        n = int(input("Enter number of nucleotides:"))
        t = int(input("Enter number of sequences:"))
        nucleotides = ['A', 'C', 'G', 'T']
        sequences = [''.join(random.choice(nucleotides) for _ in range(n)) for _ in range(t)]
        multi_aligned_seq = np.array([list(seq) for seq in sequences])
        return multi_aligned_seq, t, n


    '''Get frequency matrix'''


    def get_freq_matrix(multi_aligned_Seq, n, t):
        A_count = 0
        T_count = 0
        C_count = 0
        G_count = 0
        A_base = []
        T_base = []
        C_base = []
        G_base = []
        for i in range(n):
            for seq in multi_aligned_Seq:
                if seq[i].upper() == 'A':
                    A_count += 1
                elif seq[i].upper() == 'T':
                    T_count += 1
                elif seq[i].upper() == 'C':
                    C_count += 1
                elif seq[i].upper() == 'G':
                    G_count += 1
                else:
                    return None

            A_base.append(A_count / t)
            T_base.append(T_count / t)
            C_base.append(C_count / t)
            G_base.append(G_count / t)
            A_count = 0
            T_count = 0
            C_count = 0
            G_count = 0

        freq_mat = np.array([A_base, T_base, C_base, G_base])
        freq_mat = np.around(freq_mat, decimals=2)
        print('Frequency Matrix:\n{}'.format(freq_mat))
        print('##############################################################################')
        total_bases = t * n
        overall = []
        for seq in freq_mat:
            total = 0
            for val in seq:
                total += (val * t)
            total = round(total / total_bases, 2)
            overall.append(total)
        total_overall = sum(overall)
        # All the steps are right, I tried it on paper, I did that as the total_overall sometimes become 1.01 or 1.00001 as I rounded all the numbers
        print(total_overall)
        if ((total_overall == 1 and total_overall < 1.01) or
                (total_overall < 1 and total_overall > 0.97)):
            return freq_mat, overall
        return None


    '''Apply normalization and log, then pssm'''


    def PSSM(multi_aligned_Seq, n, t):
        freq_matrix, overall = get_freq_matrix(multi_aligned_Seq, n, t)
        normalized_freq = np.array([freq / overall[i] for i, freq in enumerate(freq_matrix)])
        normalized_freq = np.around(normalized_freq, decimals=2)
        print("Normalized frequency matrix:\n {}".format(normalized_freq))
        print('##############################################################################')
        mask = normalized_freq != 0
        pssm = np.zeros_like(normalized_freq)
        pssm[mask] = np.log2(normalized_freq[mask])
        pssm = np.around(pssm, decimals=2)
        return pssm


    '''Check if the input sequence in this family of multiple sequences or not'''


    def check_family(seq, pssm, multi_aligned_seq_matrix):
        seq = list(seq)
        weight = 0
        for i, c in enumerate(seq):
            c = c.upper()
            if c == 'A':
                weight += pssm[0][i]
            elif c == 'T':
                weight += pssm[1][i]
            elif c == 'C':
                weight += pssm[2][i]
            elif c == 'G':
                weight += pssm[3][i]
            else:
                return None
        weight = round(weight, 2)

        if weight > 0:
            print('The new sequence can be classified as a member of the sequence family.')
            seq = np.array(seq).reshape(1, -1)
            # Vertically stack the two arrays
            multi_aligned_seq_matrix = np.vstack((multi_aligned_seq_matrix, seq))
            return weight, multi_aligned_seq_matrix
        else:
            print("The sequence is non conserved sequence match.")
            return weight

    def main():
        multi_aligned_seq_matrix, t, n = get_multi_aligned_seq_mat()

        '''Print the multiple aligned sequences'''
        print('The multiple aligned Sequences matrix:\n {}'.format(multi_aligned_seq_matrix))
        print('##############################################################################')

        '''Print the PSSM'''
        pssm = PSSM(multi_aligned_seq_matrix, n, t)
        print("pssm: \n {}".format(pssm))

        '''Check if the input sequence fit to the matrix or not'''
        seq = input('Enter the sequence of length {} : '.format(n))
        try:
            if len(seq) == n:
                result = check_family(seq, pssm, multi_aligned_seq_matrix)
                if isinstance(result, tuple):
                    weight, multi_aligned_seq_matrix = result
                    print("The total match score = {}".format(weight))
                    print("the probability of the sequence fitting the matrix = {}".format(round(2 ** weight, 2)))
                    if weight > 0:
                        print('The new multiple aligned Sequences matrix:\n {}'.format(multi_aligned_seq_matrix))
                else:
                    print("The total match score = {}".format(result))
                    print("the probability of the sequence fitting the matrix = {}".format(round(2 ** result, 2)))
            else:
                raise ValueError("Sequence length must be {}".format(n))

        except ValueError as e:
            print("Invalid input:", str(e))

    main()
else:
    print("Invalid input")
