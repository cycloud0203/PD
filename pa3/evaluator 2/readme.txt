For Evaluator:
    Usage: 
    ./evaluator.sh [input.aux] [HPWL (* 1e6)] [runtime (seconds)]

    Ex: bash evaluator/evaluator.sh benchmark/ibm01/ibm01-cu85.aux 96 600
    [Eval] Your HPWL: 96 (* 1e6)
    [Eval] Your runtime: 600
    [Eval] SCORE = QUALITY_SCORE*0.8 + RUNTIME_SCORE*0.2
                = 6.9909535*0.8 + 4.1765686*0.2
                = 6.428074

    --
    Contact: Zong-Ying Cai <zycai@eda.ee.ntu.edu.tw> and Chih-Jung Hsu <alexhsu@eda.ee.ntu.edu.tw>