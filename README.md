
The entire code is contained in the file DSA_assignment_group_31.c which is in the same directory as this file.

Assumptions:
1. All input data points will have integral values for x and y (not floating point) as given in the sample data sets.
2. Any other input file used will have the same format of storing points as given in the sample input files.

Steps:
1. Open the file DSA_assignment_group_31.c
2. Scroll to the bottom to find the main function.
3. The second line in the main function will be 
        FILE* fp = fopen("./data.txt", "r");
4. Update the string that is passed as the first parameter to fopen, and specify the path of the input file you want to read the data from.
5. Two lines after that, you'll find the following loop
        for(int i=1; i<=21 && fscanf(fp, "%lf %lf\n", &x, &y)!=EOF; i++){
            Entry* e = createLeafEntry(x, y);
            Insert(e, myTree);
        }
6. In the condition of the for loop, replace 21 with the number of points you want to read from the file.
    If you wish to read all points in the file, you can remove the 'i<=21 &&' from the condition.
     (Note: printing 1 lac points takes a lot of time, took us about 10 mins)
7. Then compile and run the DSA_assignment_group_31.c
8. It will create a file named treecontent.txt in the same folder. 
9. Open treecontent.txt to see the output. 
    For internal nodes, MBRs have been printed and for leaf nodes, MBRs as well as stored 2D objects have been printed.
    (Indentations will help you identify parent-children relationships among the printed nodes easily and verify the correctness)