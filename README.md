# evolution_algorithm_parallelism_VS
PAEHN
compile envirment: Windows 10 64bit
IDE: Microsoft Visual Studio 2022
programming language: C++

option.txt: specify setting file.
paehn.txt contains the same setting as in the thesis recommended.
paehn1.txt is the variant without using recurrsive clearing when speciation.
paehn2.txt is the variant randomly producing new population.
paehn3.txt is the variant randomly updating the promising archive.
paehn4.txt is the variant which doesn't check whether the individuals in the promising archive are too crowdy and doesn't restart.

paehn_st_upper_2.txt ~ paehn_st_upper_10.txt are the paehn which use different stagnation threshold upper limit, respectively.

<font size="3"><b>OUTPUT</b></font><br />
Every time executing the paehn will produce a folder named by DateTime in the directory.
To output the coordinates and fitness value of an individual in population, promising archive and subpopulation,
please modify the kOutput parameter in "source/differential_evolution/variants/paehn.cpp" to true.
If kOutput parameter is setting to true, it will output the coordinates and fitness value of an individual in population, promising archive and subpopulation. 
Otherwise the folder will still be created, but won't contain any data.
These data can be used for graphing and interpreting.

<font size="3"><b>Graphing</b></font><br />
IDE: Matlab R2023a
Please open the "\benchmark\inc\MMOP\matlab\plots.m" by MatLab for graphing.
Modify the dirname in plots.m for graphing the specified problem.
Notice: only can graph on 1~2D functions.



編譯環境: Windows 10 64bit
IDE: Microsoft Visual Studio 2022
撰寫語言: C++

更改 option.txt 指定使用的設定檔案相對主程式路徑
paehn 為與論文參數設定相同的版本設定檔案
paehn1.txt 為不使用遞迴清除法的種化分群法
paehn2.txt 為隨機產生新主要族群
paehn3.txt 為隨機更新潛力解族群
paehn4.txt 不檢查潛力解族群是否過於擁擠且不重啟這些解個體

paehn_st_upper_2.txt ~ paehn_st_upper_10.txt 分別代表使用不同停滯閾值上限的paehn

程式每次執行會根據執行時間在目錄中產生資料夾
如欲輸出 paehn 執行過程中的主要族群、潛力解族群、子族群的座標與適應值
請修改 source/differential_evolution/variants/paehn.cpp 中的 kOutput 參數為 true，若設定 false 資料夾仍會被創建，只是裡面不會有資料
設定為 true 後便會輸出主要族群、潛力解族群、子族群的座標與適應值，可供作圖或數據判讀之用。


IDE: MatLab R2023a
作圖請使用 MatLab 開啟\benchmark\inc\MMOP\matlab\plots.m
針對想要作圖的題目修改 plots.m 中的 dirname
注意：僅能針對 1~2D的 function 作圖
