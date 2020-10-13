
if __name__ == "__main__":
    time_file = "it_0_times.txt"
    with open(f"/scratch/hartman.da/{time_file}", "r") as f:
        averages = []
        sum_time = 0
        for line in f:
            stripped = line.strip()
            if len(stripped) < 5:
                num_sims = int(stripped)
            elif len(stripped) > 22:
                averages.append(sum_time / num_sims)
                pass
            elif len(stripped) > 12:
                sum_time += float(stripped)
        print(sum(averages)/len(averages))