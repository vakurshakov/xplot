def parse(path, separator='\t'):
    results = []

    with open(path) as f:
        values = f.readline()[:-1].split(separator)
        pairs_num = len(values) // 2

        for i in range(pairs_num):
            results.append({float(values[i * 2]): float(values[i * 2 + 1])})

        for line in f.readlines():
            values = line[:-1].split(separator)

            for i in range(pairs_num):
                if values[i * 2] == '--' or values[i * 2 + 1] == '--':
                    continue

                results[i][float(values[i * 2])] = float(values[i * 2 + 1])

    return results

def parse_shared(path):
    results = []

    with open(path) as f:
        values = f.readline()[:-1].split('\t')
        pairs_num = len(values)

        for i in range(1, pairs_num):
            results.append({float(values[0]): float(values[i])})

        for line in f.readlines():
            values = line[:-1].split('\t')

            for i in range(1, pairs_num):
                results[i-1][float(values[0])] = float(values[i])

    return results