from tree import Tree
import numpy as np


def read_coef(file_name):
    """
    :param file_name:
    :return: list of lists of coefficients in a certain order
    """
    with open(file_name) as f:
        line = next(f).rstrip()  # header with version etc
        line = line.split(" ")
        coefficients = []
        for line in f:
            if line[0] == "#":
                continue
            line = line.rstrip()
            line = line.split(" ")
            coefficients.append([float(value) for value in line])
        return coefficients


def read_migration(file_name):
    """
    :param file_name:
    :return: migration rates matrix (in format: list of lists)
    """
    with open(file_name) as f:
        line = next(f).rstrip()  # header with version etc
        line = line.split(" ")
        migration_rates = []
        for line in f:
            if line[0] == "#":
                continue
            line = line.rstrip()
            line = line.split(" ")
            migration_rates.append([float(value) for value in line])
        for i in range(len(migration_rates)):
            migration_rates[i][i] = 0.0
        return(migration_rates)


def read_tree_Newick(file_name):
    """
    :param file_name:
    :return: string with tree in Newick format
    """
    with open(file_name) as f:
        line = next(f).rstrip()  # header with version etc
        line = line.split(" ")
        tree = []
        for line in f:
            if line[0] == "#":
                continue
            return line


def system_of_DE_for_lines(data: Tree, p: np.ndarray) -> np.ndarray:
    """
    :param data:
    :param p: p means P(L_i = l_i | T)
    :return: (result) derivative of probability function for each line in any population,
             1-D array with length = m * n
    """
    result = []
    value_of_cur_fun = 0
    # по каждой популяции
    for pop in range(data.number_of_populations):
        # по всем образцам
        for i in range(data.cur_samples_amount):
            # *** обращение к нужному p через индекс ( pop * data.cur_samples_amount + i )
            ################################################
            # первое(1-я часть семмы) слагаемое
            for a in range(data.number_of_populations):
                value_of_cur_fun += data.migration_probability[a][pop] * p[a * data.cur_samples_amount + i]
            ################################################
            # первое(2-я часть суммы) слагаемое
            for a in range(data.number_of_populations):
                value_of_cur_fun -= data.migration_probability[pop][a] * p[pop * data.cur_samples_amount + i]
            ################################################
            # второе слагаемое
            first_multiplier = 0
            second_multiplier = 0
            sum_of_mult = 0
            for a in range(data.number_of_populations):
                first_multiplier = data.coalescence_probability[a] * p[a * data.cur_samples_amount + i]
                second_multiplier = 0
                for k in range(data.cur_samples_amount):
                    if k != i:
                        second_multiplier += p[a * data.cur_samples_amount + k]
                sum_of_mult += first_multiplier * second_multiplier
            value_of_cur_fun += p[pop * data.cur_samples_amount + i] * sum_of_mult
            ################################################
            # третье слагаемое
            tmp_sum = 0
            for k in range(data.cur_samples_amount):
                if k != i:
                    tmp_sum += p[pop * data.cur_samples_amount + k]
            value_of_cur_fun -= p[pop * data.cur_samples_amount + i] * data.coalescence_probability[pop] * tmp_sum
            ################################################
            result.append(value_of_cur_fun)
            value_of_cur_fun = 0

    return result


def system_of_DE_for_lines(data: Tree, p: np.ndarray) -> np.ndarray:
    """
    :param data: 
    :param p: p means P(L_i = l_i | T)
    :return: (result) derivative of probability function for each line in any population,
             1-D array with length = m * n
    """
    result = []
    value_of_cur_fun = 0
    # по каждой популяции
    for pop in range(data.number_of_populations):
        # по всем образцам
        for i in range(data.cur_samples_amount):
            # *** обращение к нужному p через индекс ( pop * data.cur_samples_amount + i )
            ################################################
            # первое(1-я часть семмы) слагаемое
            for a in range(data.number_of_populations):
                value_of_cur_fun += data.migration_probability[a][pop] * p[a * data.cur_samples_amount + i]
            ################################################
            # первое(2-я часть суммы) слагаемое
            for a in range(data.number_of_populations):
                value_of_cur_fun -= data.migration_probability[pop][a] * p[pop * data.cur_samples_amount + i]
            ################################################
            # второе слагаемое
            first_multiplier = 0
            second_multiplier = 0
            sum_of_mult = 0
            for a in range(data.number_of_populations):
                first_multiplier = data.coalescence_probability[a] * p[a * data.cur_samples_amount + i]
                second_multiplier = 0
                for k in range(data.cur_samples_amount):
                    if k != i:
                        second_multiplier += p[a * data.cur_samples_amount + k]
                sum_of_mult += first_multiplier * second_multiplier
            value_of_cur_fun += p[pop * data.cur_samples_amount + i] * sum_of_mult
            ################################################
            # третье слагаемое
            tmp_sum = 0
            for k in range(data.cur_samples_amount):
                if k != i:
                    tmp_sum += p[pop * data.cur_samples_amount + k]
            value_of_cur_fun -= p[pop * data.cur_samples_amount + i] * data.coalescence_probability[pop] * tmp_sum
            ################################################
            result.append(value_of_cur_fun)
            value_of_cur_fun = 0

    return result


def create_initial(data: Tree,
                   period: int,
                   previous_states: np.ndarray = None,
                   lineage: np.ndarray = None) -> np.ndarray:
    """
    :param data:
    :param period:
    :param previous_states: probabilities before coalescence
    :param lineage: contains two samples that participated in coalescence
    :return: result: vector of initial states for current "period"
    """
    # *** обращение к нужному p через индекс ( pop * data.cur_samples_amount + i )
    if period == 0:
        return data.get_initial_states()
    else:
        lineage = np.sort(lineage)
        result = []
        # cur_samples_sum = 0  ??нужно??

        # посчитаю все условные вероятности для lineage[0] и всех 'pop'
        conditional_prob = []
        for pop in range(data.number_of_populations):
            conditional_prob.append(previous_states[pop * (data.cur_samples_amount + 1) + lineage[0]] *
                                    previous_states[pop * (data.cur_samples_amount + 1) + lineage[1]] *
                                    data.coalescence_probability[pop])
        sum_conditional_prob = np.sum(conditional_prob)  # ??? нужно ли переводить в массив ???
        # по каждой популяции
        for pop in range(data.number_of_populations):
            # по каждому образцу
            for i in range(data.cur_samples_amount + 1):
                # если образец участвовал в коалесценции и у него меньшее 'id', то есть он остался
                if i == lineage[0]:
                    # !!!!!!!!!!!!!! можно использовать уже посчитанные выше
                    result.append(previous_states[pop * (data.cur_samples_amount + 1) + i] *
                                  previous_states[pop * (data.cur_samples_amount + 1) + lineage[1]] *
                                  data.coalescence_probability[pop])
                # если образец участвовал в коалесценции и у него большее 'id', то есть он не остался
                elif i == lineage[1]:  # добавила для лучшего понимания кода, но вообще этот if можно убрать
                    continue
                # если образец не участвовал в коалесценции
                elif i != lineage[0] and i != lineage[1]:
                    result.append(previous_states[pop * (data.cur_samples_amount + 1) + i] * sum_conditional_prob)

            # cur_samples_sum += tree.number_of_samples[pop]  ??нужно??
        return result


def DE_for_tree(data: Tree,
                p: float,
                lines_prob: np.ndarray,
                time: int,
                period_start_time: int) -> np.ndarray:
    """
    :param data:
    :param p: P_{t}(T)
    :param lines_prob: 2-D array with shape (data.cur_samples_amount x duration of current period)
                       the probability of every lines being in any state at each point of time in the current period
    :param time: specifies which column of lines_prob to use
    :param period_start_time: specifies the start time of the period to subtract from 't'
    :return: derivative of the tree probability function
    """
    sum_of_mult = 0
    cur_lines_prob = lines_prob[:, int(time - period_start_time)]
    # по каждой популяции
    for pop in range(data.number_of_populations):
        # по всем образцам
        for i in range(data.cur_samples_amount):
            # по всем образцам
            for j in range(data.cur_samples_amount):
                sum_of_mult += (0.5 * data.coalescence_probability[pop] *
                                cur_lines_prob[pop * data.cur_samples_amount + i] *
                                cur_lines_prob[pop * data.cur_samples_amount + j])

    return (-1) * p * sum_of_mult


def create_init_state_for_tree(data: Tree,
                               period: int,
                               p_tree_before: float = None,
                               lines_prob: np.ndarray = None,
                               lineage: np.ndarray = None) -> float:
    """
    :param data:
    :param period:
    :param lines_prob: 2-D array with shape (data.cur_samples_amount x 1)
                       the probability of every lines being in any state at time before the coalescence
    :param lineage: contains two samples that participated in coalescence
    :return: the probability in form of 1-D array with shape (1 x 1)
    """
    if period == 0:
        # как-то ставится задача Коши
        return np.array([1])  # !!!!!!!!!!!!!!! пока что не знаю, как ставить задачу Коши на первой итерации !!!!!!!!!!!!!!!
    else:
        # вероятность дерева после коалесценции P_{t}^{after}
        sum_of_mult = 0
        # lines_prob = lines_prob.transpose()
        for pop in range(data.number_of_populations):
            sum_of_mult += (data.coalescence_probability[pop] *
                            lines_prob[pop * data.cur_samples_amount + lineage[0]] *
                            lines_prob[pop * data.cur_samples_amount + lineage[1]])
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!надо ли делать (data.cur_sample_amount + 1) ???????????????

        return np.array([p_tree_before * sum_of_mult])
