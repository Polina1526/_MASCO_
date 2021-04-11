import numpy as np
from scipy import special   # для биномиальных коэффициентов
import struct


class Tree:

    def __init__(self, file_name: str):
        # Извлечение данных из файла и создание объекта
        # Данные в файле должны быть в таком порядке: N, k, n, t, m, q, Q
        with open(file_name, 'r', encoding='utf8') as file:
            tmp_data = []
            for line in file:
                n1 = line.find('=')
                tmp_data.append(line[n1 + 2:-1:])
            N, k, n, t, m, q, Q, tree = tmp_data
            # для отладки
            # print(N, k, n, t, m, q, Q)
        self.__N = int(N)      # эталонный размер популяции
        self.__number_of_populations = int(k)  # количество популяций (m)
        self.__number_of_samples = np.array(list(map(int, n.split(' '))))
        self.__samples_amount = np.sum(self.__number_of_samples)  #всего образцов дано (n)
        self.__T = int(t)      # время слияния всех популяций
        self.__migration_probability = np.array(list(map(float, m.split(' '))))
        self.__migration_probability.resize(int(k), int(k))
        self.__coalescence_probability = np.array(list(map(float, q.split(' '))))
        self.__Q = float(Q)
        self.__tree_newick = tree  # норм, что строка?

        self.__cur_samples_amount = self.__samples_amount

    def show(self):
        print(self.__number_of_populations, self.__T, '\n')
        print(self.__number_of_samples, self.__samples_amount, '\n')
        print(self.__migration_probability, '\n\n', self.__coalescence_probability, '\n')
        print(self.__tree_newick, '\n')

    @property
    def original_size(self):
        return self.__N

    @property
    def number_of_populations(self):
        return self.__number_of_populations

    @property
    def number_of_samples(self):
        return self.__number_of_samples

    @property
    def samples_amount(self):
        return self.__samples_amount

    @property
    def T(self):
        return self.__T

    @property
    def migration_probability(self):
        return self.__migration_probability

    @property
    def coalescence_probability(self):
        return self.__coalescence_probability

    @property
    def Q(self):
        return self.__Q

    @property
    def tree_newick(self):
        return self.__tree_newick

    @property
    def cur_samples_amount(self):
        return self.__cur_samples_amount

    @cur_samples_amount.setter
    def cur_samples_amount(self, cur_samples_amount):
        if cur_samples_amount >= 1:
            self.__cur_samples_amount = cur_samples_amount
        else:
            raise ValueError("the minimum number of samples has been reached")

    def get_initial_states(self):
        """
        :return result: 1-D array of initial states
        """
        result = []
        # вспомогательная переменная для проверки принадлежности текущего образца к текущей популяции
        current_samples_sum = 0
        # по каждой популяции
        for i in range(self.__number_of_populations):
            # по каждому образцу
            for j in range(self.__samples_amount):
                # если образец в начальный момент принадлежит этой популяции
                if(current_samples_sum <= j < current_samples_sum + self.__number_of_samples[i]):
                    result.append(1)
                else:
                    result.append(0)
            current_samples_sum += self.__number_of_samples[i]
        return result
