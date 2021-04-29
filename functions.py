from tree import Tree
import numpy as np
from struct import unpack


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
            # первое(1-я часть суммы) слагаемое
            for a in range(data.number_of_populations):
                value_of_cur_fun += data.migration_probability[a][pop] * p[a * data.cur_samples_amount + i]
            ################################################
            # первое(2-я часть суммы) слагаемое
            for a in range(data.number_of_populations):
                value_of_cur_fun -= data.migration_probability[pop][a] * p[pop * data.cur_samples_amount + i]
            ################################################
            # второе слагаемое
            sum_of_mult = 0
            for a in range(data.number_of_populations):
                tree_prob_mult = 0
                for j in range(data.cur_samples_amount):
                    if j != i:
                        for k in range(data.cur_samples_amount):
                            if (k != i) and (k != j):
                                # вероятность дерева по k-ым линиям
                                tree_sum_k = 0
                                # вероятность дерева по j-ым линиям
                                tree_sum_j = 0
                                for l in range(data.number_of_populations):
                                    tree_sum_k += p[l * data.cur_samples_amount + k]
                                    tree_sum_j += p[l * data.cur_samples_amount + j]
                                # сумма произведений всех условных вероятностей
                                tree_prob_mult += ((p[pop * data.cur_samples_amount + k] / tree_sum_k) *
                                                   (p[pop * data.cur_samples_amount + j] / tree_sum_j))
                sum_of_mult += data.coalescence_probability[a] * tree_prob_mult
            value_of_cur_fun -= 0.5 * p[pop * data.cur_samples_amount + i] * sum_of_mult
            """
            space
            """
            ################################################
            # третье слагаемое
            cond_prob = 0
            for k in range(data.cur_samples_amount):
                if k != i:
                    # расчёт вероятности дерева, для подсчёта условной вероятности
                    tree_prob = 0
                    for l in range(data.number_of_populations):
                        tree_prob += p[l * data.cur_samples_amount + k]
                    # условная вероятность
                    cond_prob += p[pop * data.cur_samples_amount + k] / tree_prob
            value_of_cur_fun -= p[pop * data.cur_samples_amount + i] * data.coalescence_probability[pop] * cond_prob
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
        # посчитаю все вероятности для lineage[0] и всех 'pop'
        line_prob_after = []
        # вероятность дерева по i-ым(lineage[0]) линиям
        tree_sum_i = 0
        # вероятность дерева по j-ым(lineage[1]) линиям
        tree_sum_j = 0
        for l in range(data.number_of_populations):
            tree_sum_i += previous_states[l * (data.cur_samples_amount + 1) + lineage[0]]
            tree_sum_j += previous_states[l * (data.cur_samples_amount + 1) + lineage[1]]
        for pop in range(data.number_of_populations):
            line_prob_after.append(previous_states[pop * (data.cur_samples_amount + 1) + lineage[0]] *
                                    (previous_states[pop * (data.cur_samples_amount + 1) + lineage[1]] / tree_sum_j) *
                                    data.coalescence_probability[pop])
        sum_line_prob_after = np.sum(line_prob_after)
        # по каждой популяции
        for pop in range(data.number_of_populations):
            # по каждому образцу
            for i in range(data.cur_samples_amount + 1):
                # если образец участвовал в коалесценции и у него меньшее 'id', то есть он остался
                if i == lineage[0]:
                    # !!!!!!!!!!!!!! можно использовать уже посчитанные выше
                    result.append(previous_states[pop * (data.cur_samples_amount + 1) + i] *
                                  (previous_states[pop * (data.cur_samples_amount + 1) + lineage[1]] / tree_sum_j) *
                                  data.coalescence_probability[pop])
                # если образец участвовал в коалесценции и у него большее 'id', то есть он не остался
                elif i == lineage[1]:  # добавила для лучшего понимания кода, но вообще этот if можно убрать
                    continue
                # если образец не участвовал в коалесценции
                elif i != lineage[0] and i != lineage[1]:
                    # вероятность дерева по k-ым линиям
                    tree_sum_k = 0
                    for l in range(data.number_of_populations):
                        tree_sum_k += previous_states[l * (data.cur_samples_amount + 1) + i]
                    result.append((previous_states[pop * (data.cur_samples_amount + 1) + i] / tree_sum_k) *
                                  sum_line_prob_after)

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
                if j != i:
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
    :param p_tree_before:
    :param lines_prob: 2-D array with shape (data.cur_samples_amount x 1)
                       the probability of every lines being in any state at time before the coalescence
    :param lineage: contains two samples that participated in coalescence
    :return: the probability in form of 1-D array with shape (1 x 1)
    """
    if period == 0:
        # в самом начале вероятность дерева равна 1
        return np.array([1])
    else:
        # вероятность дерева после коалесценции P_{t}^{after}
        sum_of_mult = 0
        for pop in range(data.number_of_populations):
            sum_of_mult += (data.coalescence_probability[pop] *
                            lines_prob[pop * data.cur_samples_amount + lineage[0]] *
                            lines_prob[pop * data.cur_samples_amount + lineage[1]])

        return np.array([p_tree_before * sum_of_mult])


def parse_tree_history(data: Tree) -> tuple:
    """
    :param data:
    :return: tuple of the lists with limits of integration and lines that coalesced (limits_list, lineage_list)
    """
    limits_list = []
    lineage_list = []
    shift = np.zeros(data.samples_amount, dtype=int)

    out_file_bin = open("tree_history.bin", "rb")

    # нахождение количества событий
    out_file_bin.seek(0, 2)
    file_size = out_file_bin.tell()
    event_size = 4 * 5
    event_amount = int((file_size - 4) / (4 * 5))  # так как в начале записано TMRCA

    # переход в начало файла
    out_file_bin.seek(4, 0)
    # по каждому событию
    for i in range(event_amount):
        event = (unpack('<5i', out_file_bin.read(event_size)))
        # пропускаются события миграции
        if event[0] == 1:
            continue
        # обработка события коалесценции
        else:
            ev, pop_id, offspr1, offspr2, time = event
            if len(limits_list) == 0:
                limits_list.append(np.array([0, time]))
            else:
                prev_time = limits_list[-1][1]
                limits_list.append(np.array([prev_time, time]))

        lineage_list.append(np.array([offspr1 - shift[offspr1], offspr2 - shift[offspr2]]))
        # устанановка сдвигов индексов, так как образцы перенумеровываются после коалесценции
        for j in range(offspr2 + 1, len(shift)):
            shift[j] += 1

    out_file_bin.close()

    return limits_list, lineage_list
