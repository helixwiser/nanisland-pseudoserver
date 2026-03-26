#!/usr/bin/env python3




import sys
import csv
from collections import defaultdict




def main(input_file, output_file):
    # 存储每个 organism 的指标列表
    stats = defaultdict(list)




    # 读取输入文件
    with open(input_file, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            organism = row.get("Organism", "").strip()
            if not organism:
                continue




            try:
                degree = float(row["Degree"])
                deg_centrality = float(row["Degree_Centrality"])
                betweenness = float(row["Betweenness"])
                closeness = float(row["Closeness"])
                pagerank = float(row["PageRank"])
            except (ValueError, TypeError):
                continue  # 跳过无效数据




            stats[organism].append({
                "degree": degree,
                "deg_centrality": deg_centrality,
                "betweenness": betweenness,
                "closeness": closeness,
                "pagerank": pagerank
            })




    # 计算平均值并写入 CSV
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        # 写表头
        writer.writerow([
            "Organism",
            "Degree",
            "Degree_Centrality",
            "Betweenness",
            "Closeness",
            "PageRank"
        ])




        # 按 organism 排序输出
        for organism in sorted(stats.keys()):
            values = stats[organism]
            n = len(values)




            avg_degree = sum(v["degree"] for v in values) / n
            avg_centrality = sum(v["deg_centrality"] for v in values) / n
            avg_betweenness = sum(v["betweenness"] for v in values) / n
            avg_closeness = sum(v["closeness"] for v in values) / n
            avg_pagerank = sum(v["pagerank"] for v in values) / n




            writer.writerow([
                organism,
                avg_degree,
                avg_centrality,
                avg_betweenness,
                avg_closeness,
                avg_pagerank
            ])




    # 打印汇总信息到屏幕
    print("Organism_Average_Statistics")
    print("Organism,Degree,Degree_Centrality,Betweenness,Closeness,PageRank")




    for organism in sorted(stats.keys()):
        values = stats[organism]
        n = len(values)
        avg_degree = sum(v["degree"] for v in values) / n
        avg_centrality = sum(v["deg_centrality"] for v in values) / n
        avg_betweenness = sum(v["betweenness"] for v in values) / n
        avg_closeness = sum(v["closeness"] for v in values) / n
        avg_pagerank = sum(v["pagerank"] for v in values) / n




        print(organism + "," +
              str(avg_degree) + "," +
              str(avg_centrality) + "," +
              str(avg_betweenness) + "," +
              str(avg_closeness) + "," +
              str(avg_pagerank))




    print("\nSummary")
    print("Total_Organisms," + str(len(stats)))
    total_nodes = sum(len(vals) for vals in stats.values())
    print("Total_Nodes," + str(total_nodes))




    print("\n✅ Statistics saved to: " + output_file)




if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python analyze_organism_stats.py <input.tsv> <output.csv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])