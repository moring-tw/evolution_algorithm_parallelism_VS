# evolution_algorithm_parallelism_VS
PAEHN<br />
compile envirment: Windows 10 64bit<br />
IDE: Microsoft Visual Studio 2022<br />
programming language: C++<br />

option.txt: specify setting file.<br />
paehn.txt contains the same setting as in the thesis recommended.<br />
paehn1.txt is the variant without using recurrsive clearing when speciation.<br />
paehn2.txt is the variant randomly producing new population.<br />
paehn3.txt is the variant randomly updating the promising archive.<br />
paehn4.txt is the variant which doesn't check whether the individuals in the promising archive are too crowdy and doesn't restart.<br />

paehn_st_upper_2.txt ~ paehn_st_upper_10.txt are the paehn which use different stagnation threshold upper limit, respectively.<br />

<font size="3"><b>OUTPUT</b></font><br />
Every time executing the paehn will produce a folder named by DateTime in the directory.<br />
To output the coordinates and fitness value of an individual in population, promising archive and subpopulation,<br />
please modify the kOutput parameter in "source/differential_evolution/variants/paehn.cpp" to true.<br />
If kOutput parameter is setting to true, it will output the coordinates and fitness value of an individual in population, promising archive and subpopulation. <br />
Otherwise the folder will still be created, but won't contain any data.<br />
These data can be used for graphing and interpreting.<br />

<font size="3"><b>Graphing</b></font><br />
IDE: Matlab R2023a<br />
Please open the "\benchmark\inc\MMOP\matlab\plots.m" by MatLab for graphing.<br />
Modify the dirname in plots.m for graphing the specified problem.<br />
Notice: only can graph on 1~2D functions.<br />



編譯環境: Windows 10 64bit<br />
IDE: Microsoft Visual Studio 2022<br />
撰寫語言: C++<br />

更改 option.txt 指定使用的設定檔案相對主程式路徑<br />
paehn 為與論文參數設定相同的版本設定檔案<br />
paehn1.txt 為不使用遞迴清除法的種化分群法<br />
paehn2.txt 為隨機產生新主要族群<br />
paehn3.txt 為隨機更新潛力解族群<br />
paehn4.txt 不檢查潛力解族群是否過於擁擠且不重啟這些解個體<br />

paehn_st_upper_2.txt ~ paehn_st_upper_10.txt 分別代表使用不同停滯閾值上限的paehn<br />

程式每次執行會根據執行時間在目錄中產生資料夾<br />
如欲輸出 paehn 執行過程中的主要族群、潛力解族群、子族群的座標與適應值<br />
請修改 source/differential_evolution/variants/paehn.cpp 中的 kOutput 參數為 true，若設定 false 資料夾仍會被創建，只是裡面不會有資料<br />
設定為 true 後便會輸出主要族群、潛力解族群、子族群的座標與適應值，可供作圖或數據判讀之用。<br />


IDE: MatLab R2023a<br />
作圖請使用 MatLab 開啟\benchmark\inc\MMOP\matlab\plots.m<br />
針對想要作圖的題目修改 plots.m 中的 dirname<br />
注意：僅能針對 1~2D的 function 作圖<br />
