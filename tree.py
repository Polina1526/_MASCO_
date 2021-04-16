import numpy as np
import functions as fun
#from functions import read_coef, read_migration, read_tree_Newick
from scipy import special   # для биномиальных коэффициентов
import struct


class Tree:
    def __init__(self, coef_file: str, migr_file: str, tree_file: str):
        # Извлечение данных из файла и создание объекта
        # Данные в файле должны быть в таком порядке: N, k, n, t, m, q, Q
        """with open(file_name, 'r', encoding='utf8') as file:
            tmp_data = []
            for line in file:
                n1 = line.find('=')
                tmp_data.append(line[n1 + 2:-1:])
            N, k, n, t, m, q, Q, tree = tmp_data
            # для отладки
            # print(N, k, n, t, m, q, Q)"""
        N, k, n, t, q, Q = fun.read_coef(coef_file)
        self.__N = int(N[0])      # эталонный размер популяции
        self.__number_of_populations = int(k[0])  # количество популяций (m)
        self.__number_of_samples = np.array(n)  # (список) колличество образцов в каждой популяции
        self.__samples_amount = int(np.sum(self.__number_of_samples))  # всего образцов (n)
        self.__T = int(t[0])      # время слияния всех популяций
        self.__migration_probability = np.array(fun.read_migration(migr_file))  # матрица с коэфициентами миграции
        self.__coalescence_probability = np.array(q)  # вектор с коэффициентами коалесценции
        self.__Q = float(Q[0])
        self.__tree_newick = fun.read_tree_Newick(tree_file)  # дерево в формате Newick норм, что строка?

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
