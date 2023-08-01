import math
import re
import statistics

from sklearn.linear_model import LinearRegression

NA_VALUE = "NA"
ID, SEQ, CY3, CY3_FLAGS, ALEXA, ALEXA_FLAGS = range(1, 7)


def _get_margins(array_type, func_name):
    if array_type == "4x44k":
        return 20, 245, 150, 416
    if array_type == "8x15k":
        return 10, 86, 154, 250
    raise ValueError(f"Array type {array_type} is not supported for {func_name}")


def _get_top_row_col(array_type, func_name):
    if array_type == "4x44k":
        return 170, 266
    if array_type == "8x15k":
        return 164, 96
    raise ValueError(f"Array type {array_type} is not supported for {func_name}")


def median_neighborhood(data_matrix_Aref, radius, fluor, array_type, chip_median, last_cell):
    """
    For each spot on a 4x44K or 8x15K array, takes the neighborhood
        of spots as the block of radius R (up, down, left,
        right by R rows and columns), and computes the median
        intensity for the block.
    4x44k array has 45220 spots: 170 rows x 266 columns.
    8x15K array has 16032 spots: 167 rows x 96 columns.
    Adjusts for rounded corners and margins (11/20/06).  If
        the spot is too close to the corner, go in further
        (diagonally) so the block size is constant.  If the
        spot is too close to the margin, go in further
        (laterally) so the block size is constant.

    M. Berger
    """
    next_cell = last_cell + 1

    if fluor == "Cy3":
        adj_bsi, flags = CY3, CY3_FLAGS
    elif fluor == "Alexa488":
        adj_bsi, flags = last_cell, ALEXA_FLAGS
    else:
        raise ValueError("Incorrect usage of median_neighborhood_4x44k")

    upper_row, upper_col = _get_top_row_col(array_type, "median_neighborhood")
    upper_left, upper_right, lower_left, lower_right = _get_margins(array_type, "median_neighborhood")

    for row in range(1, upper_row+1):
        for col in range(1, upper_col+1):
            cell = data_matrix_Aref[col][row]
            # Upper Left Corner: Empty spots defined by the line "Row + Column < upper_left".
            if col + row < upper_left + 2 * radius:
                counter = math.ceil((upper_left + 2 * radius - row - col) / 2)
                adj_col, adj_row = col + counter, row + counter
            # Upper Right Corner: Empty spots defined by the line "Column - Row > upper_right".
            elif col - row > upper_right - 2 * radius:
                counter = math.ceil((-upper_right + 2 * radius + col - row) / 2)
                adj_col, adj_row = col - counter, row + counter
            # Lower Left Corner: Empty spots defined by the line "Row - Column > lower_left".
            elif row - col > lower_left - 2 * radius:
                counter = math.ceil((-lower_left + 2 * radius + row - col) / 2)
                adj_row, adj_col = row - counter, col + counter
            # Lower Right Corner: Empty spots defined by the line "Row + Column > lower_right".
            elif row + col > lower_right - 2 * radius:
                counter = math.ceil((-lower_right + 2 * radius + row + col) / 2)
                adj_row, adj_col = row - counter, col - counter
            # All other spots
            else:
                adj_col, adj_row = col, row

            # Row, Column Boundaries
            def boundaries(adj, upper, radius):
                if radius < adj <= upper - radius:
                    return adj - radius, adj + radius
                if adj <= radius:
                    return 1, 1 + 2 * radius
                if adj > upper - radius:
                    return upper - 2 * radius, upper

            # Row, Column Boundaries
            top_row, bottom_row = boundaries(adj_row, upper_row, radius)
            left_col, right_col = boundaries(adj_col, upper_col, radius)

            # Make array of elements in block
            block_spots = []
            for block_col in range(left_col, right_col+1):
                for block_row in range(top_row, bottom_row+1):
                    block_cell = data_matrix_Aref[block_col][block_row]
                    if (block_cell[adj_bsi] != NA_VALUE and block_cell[flags] > -100 and
                        (block_cell[ID].startswith("dBr") or block_cell[ID].startswith("Ctrl"))):
                        block_spots.append(block_cell[adj_bsi])

            block_size = (2 * radius + 1) * (2 * radius + 1)

            if len(block_spots) <= block_size / 2:
                cell.append(chip_median)
            else:
                cell.append(median(block_spots))

            if cell[flags] > -100 and fluor == "Cy3" and \
                    (cell[ID].startswith("dBr") or cell[ID].startswith("Ctrl")):
                cell.append(cell[adj_bsi] / cell[next_cell] * chip_median)
            elif cell[flags] > -100 and fluor == "Alexa488" and cell[last_cell] != NA_VALUE:
                cell.append(cell[adj_bsi] / cell[next_cell] * chip_median)
            else:
                data_matrix_Aref[col][row].append(NA_VALUE)

    return next_cell + 1


def cy3_regression(data_matrix, regression_order, array_type, output, last_cell):
    """
    Computes the regression coefficients for Cy3 using the (unflagged)
        combinatorial "dBr" spots.  Appends expected Cy3 to data matrix.
        Calculates observed/expected Cy3 for all custom-designed combinatorial
        and control spots.
            -- adapted from cpan.org, A. Philippakis, M. Berger
    """
    top_row, top_col = _get_top_row_col(array_type, 'cy3_regression')
    combinatorial_length = regression_order + 36  #MFB changed 2/10/07

    x_component_names = ["A" + convert_to_letters(i, regression_order) for i in range(4 ** regression_order)]

    X, y = [], []
    signal_cell = last_cell if last_cell > 6 else 3

    for row in range(1, top_row + 1):
        for col in range(1, top_col + 1):
            cell = data_matrix[col][row]
            if "dBr" in cell[ID] and cell[last_cell] != NA_VALUE and cell[CY3_FLAGS] > -100:
                sequence = cell[SEQ][:combinatorial_length]
                X.append([sum(1 for _ in re.finditer(f'(?=({key}))', sequence, re.IGNORECASE))
                     for key in x_component_names])
                y.append(cell[signal_cell])

    # Perform linear regression
    reg = LinearRegression()
    reg.fit(X, y)
    theta = [reg.intercept_] + list(reg.coef_)

    with open(output, "w") as f:
        f.write("*****************************************************\n"
                "Regression 'components'\n"
                "*****************************************************\n")
        f.writelines([f"Theta[{i}='{mer}']=\t{weight}\n"
                      for i, (mer, weight) in enumerate(zip(['intercept'] + x_component_names, theta))])

        r_squared = round(reg.score(X, y), 3)
        f.write(f"R^2 = {r_squared}\n"
                "*****************************************************\n")

    for row in range(1, top_row + 1):
        for col in range(1, top_col + 1):
            cell = data_matrix[col][row]
            if ("dBr" in cell[ID] or "Ctrl" in cell[ID]) and cell[last_cell] != NA_VALUE and cell[CY3_FLAGS] > -100:
                array_sequence = cell[SEQ][:combinatorial_length]
                x = [1] + [
                    sum(1 for _ in re.finditer(f'(?=({key}))', array_sequence, re.IGNORECASE))
                    for key in x_component_names
                ]

                predicted_value = sum(mer * weight for (mer, weight) in zip(x, theta))
                cell.append(predicted_value)
                cell.append(cell[signal_cell] / predicted_value)
            else:
                cell.extend([NA_VALUE, NA_VALUE])

    return last_cell + 2


def alexa488_normalize(data_matrix_Aref, alexa488_list, array_type, last_cell):
    """
    Normalizes Alexa488 signal by (Observed/Expected) Cy3 and stores
        array of all Alexa488 intensities for custom-designed combinatorial
        and control probes.
    """
    top_row, top_col = _get_top_row_col(array_type, 'alexa488_normalize')

    for row in range(1, top_row+1):
        for col in range(1, top_col+1):
            cell = data_matrix_Aref[col][row]
            if cell[last_cell] != NA_VALUE and 0.5 < cell[last_cell] < 2 and cell[ALEXA_FLAGS] > -100:
                cell.append(cell[ALEXA] / cell[last_cell])
                alexa488_list.append(cell[-1])
            else:
                cell.append(NA_VALUE)

    return last_cell + 1


def median(numbers):
    """
    Returns the median of an array of numbers
    """
    return statistics.median(numbers)


def convert_to_letters(num, k):
    """
    Takes a number and a "k" and returns the ACGT-equivalent of that number for length k
        -- written by A. Philippakis
    """
    return ''.join('ACGT'[(num >> i) & 3] for i in range((k - 1) * 2, -1, -2))
