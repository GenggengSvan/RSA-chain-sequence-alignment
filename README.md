# RSA-chain-sequence-alignment
本题抽象为输入两个病毒的 RNA 个数和每一条 RNA 的碱基序列，进行相似性比较。如：有两个不同病毒 A 和 B，A 病毒有 3 条 RNA 序列，分别为：①ACUUCG ②GUGAUAGACAC③AUGGAGA；B 病毒有 5 条 RNA 序列，分别为：①ACUAG ②GUGAUAGACAC ③AUCGUGUG④CCGUGGA ⑤CUAUGUG。在进行相似性比对时，需要对 A 病毒和 B 病毒的多条 RNA 进行匹配，如 A 病毒的 RNA①和 B 病毒的 RNA①进行匹配，这两条 RNA 匹配的得分由序列比对这两条 RNA 的碱基序列的方式得出。需要寻找一个最优的匹配方式，使匹配的总得分最高，由这个得分代表两个病毒的相似度。
使用序列比对算法进行RSA链比对，性能一般，仅供参考。
