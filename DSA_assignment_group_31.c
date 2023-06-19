#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#define M 4     // Maximum number of entries in a node
#define m 2     // Minimum number of entries in a node
FILE* outp;

typedef struct point Point;
typedef struct rectangle Rectangle;
typedef struct node Node;
typedef struct entry Entry;
typedef struct r_tree R_Tree;

struct point{   // Stores x and y coordinates of a 2-D point
    double x, y;
};

struct rectangle{   // Stores the top-right and bottom-left vertices of a rectangle
    Point bottom_left, top_right;
};

// Describes each entry in a node
struct entry{   
    Rectangle I;    // Minimum bounding rectangle (MBR) of the child node or of the terminal object stored
    unsigned int is_leaf_entry;     // Whether or not the entry is in a leaf node
    union{      
        /*Either an entry will store a terminal object (if it is a leaf entry) or it will store a pointer to a 
        child node (if it is non-leaf entry) but not both; therefore union is used*/
       
        Rectangle obj;  // 2-D points will be stored as rectangles with coincident bottom-left and top-right vertices
        Node* p;
    } data;
};

// Defines a node
struct node{    
    unsigned int is_leaf;   //  If it is a leaf node or not
    int occupied;           // Number of entries currently present in the node (>=m and <=M)
    Entry* E[M];            // Array of pointers to the entries of this node
};

// The R-Tree stucture
struct r_tree{  
    Node* T;    // Root node
    Entry* Meta_Entry;  // A dummy entry that has the rectangle enclosing the entries in the root node
};

// Create an entry; 'is_L' is 1 if the entry to be created is a leaf entry
Entry* createEntry(unsigned int is_L){  
    Entry* temp_entry = (Entry*) malloc(sizeof(Entry));
    temp_entry->is_leaf_entry = is_L;

    //Initialize the MBR of the new entry to a point (0,0)
    temp_entry->I.bottom_left.x =0;
    temp_entry->I.bottom_left.y =0;
    temp_entry->I.top_right.x =0;
    temp_entry->I.top_right.y =0;
    return temp_entry;
}

// Create a node; 'is_L' is 1 if the node to be created is a leaf node
Node* createNode(unsigned int is_L){    
    Node* temp_node = (Node*) malloc(sizeof(Node));
    temp_node->is_leaf = is_L;
    temp_node->occupied = 0;    //Initialize number of entries in the node to 0
    return temp_node;
}

// Create an R-Tree
R_Tree* create_R_Tree(){    
    R_Tree* temp_tree = (R_Tree*) malloc(sizeof(R_Tree));
    temp_tree->T = createNode(1);   // Create a root node right after creating the tree
    temp_tree->Meta_Entry = createEntry(0);
    return temp_tree;
}

// Checks if two rectangles overlap
unsigned int overlaps(Rectangle r1, Rectangle r2){  
    return (r1.bottom_left.x <= r2.top_right.x && r1.top_right.x >= r2.bottom_left.x &&
            r1.bottom_left.y <= r2.top_right.y && r1.top_right.y >= r2.bottom_left.y);
}

// Just prints the rectangle (which is actually a point)
void print_leaf_object(Rectangle r){    
    printf("\nLeaf Node:\n");
    printf("(%lf, %lf)", r.bottom_left.x, r.bottom_left.y);
}

// Search all objects that overlap 'search_rectangle' in the tree rooted at T
void search_R_Tree(Node* T, Rectangle search_rectangle){
    // if T is not a leaf, recursively search in all child nodes whose MBR overlaps the given rectangle    
    if(!(T->is_leaf)){          
        for(int i=0; i<T->occupied; i++){
            if(overlaps(T->E[i]->I, search_rectangle))
                search_R_Tree(T->E[i]->data.p, search_rectangle);
        }
    }
    // If T is a leaf node, print all entries in T that overlap the given rectangle
    else{   
        for(int i=0; i<T->occupied; i++){
            if(overlaps(T->E[i]->I, search_rectangle)){
                printf("\n");
                print_leaf_object(T->E[i]->data.obj);
            }
        }
    }
}

double Area(Rectangle r){   // Finds area of a rectangle
    return ((r.top_right.x)-(r.bottom_left.x))*((r.top_right.y)-(r.bottom_left.y));
}

// Checks if r1 is completely contained in (lies within) r2
int is_contained(Rectangle r1, Rectangle r2){   
    return (r1.bottom_left.x >= r2.bottom_left.x && r1.bottom_left.y >= r2.bottom_left.y &&
            r1.top_right.x <= r2.top_right.x && r1.top_right.y <= r2.top_right.y);
}

// Adds entry e into node n and increments the occupied variable of the node
void add_entry_to_node(Entry* e, Node* n){  
    n->E[(n->occupied)++] = e;
}

// Returns the most appropriate leaf node in the tree rooted at N for the insertion of entry E
Node* ChooseLeaf(Node* N, Entry* E){
    // Base case: If the root given is already a leaf node, return that node    
    if(N->is_leaf) return N;    

    int entry_index=0;
    double enlarge_reqd = 999999999, area = 999999999, temp_area_increase;

    // Assume first entry to be the appropriate one
    Entry* the_one = N->E[0];  

    for(int i=1; i<N->occupied; i++){      // Traverse through other entries of the node
        if(is_contained(E->I, N->E[i]->I)){     
            /* If there's an entry whose child node's MBR needs no enlargement to accomodate the new entry, 
            then certainly choose this and break*/
            enlarge_reqd = 0;
            area = Area(N->E[i]->I);
            entry_index = i;
            break;
        }
        else{   
            // Check if we find a better entry whose child node's MBR needs lesser enlargement than the one we assumed in 'the_one'
            Rectangle bigger_rect = N->E[i]->I;
            if(E->I.bottom_left.x < N->E[i]->I.bottom_left.x) bigger_rect.bottom_left.x = E->I.bottom_left.x;
            if(E->I.bottom_left.y < N->E[i]->I.bottom_left.y) bigger_rect.bottom_left.y = E->I.bottom_left.y;
            if(E->I.top_right.x > N->E[i]->I.top_right.x) bigger_rect.top_right.x = E->I.top_right.x;
            if(E->I.top_right.y > N->E[i]->I.top_right.y) bigger_rect.top_right.y = E->I.top_right.y;

            temp_area_increase = Area(bigger_rect) - Area(N->E[i]->I);

            if(temp_area_increase < enlarge_reqd){
                enlarge_reqd = temp_area_increase;
                area = Area(N->E[i]->I);
                entry_index = i;
            }
            // If enlargement required is same for multiple nodes, choose the one with smaller current area
            else if(temp_area_increase == enlarge_reqd && Area(N->E[i]->I) < area){ 
                area = Area(N->E[i]->I);
                entry_index = i;
            }
        }
    }
    // Assign the chosen one and recurse
    if(N->E[entry_index] && N->E[entry_index] != the_one) the_one = N->E[entry_index]; 
    return ChooseLeaf(the_one->data.p, E);
}

// Returns a rectangle that tightly encloses the two given entries
Rectangle Enclose_two(Entry* e1, Entry* e2){    
    Rectangle r = e1->I;
    if(e2->I.bottom_left.x < r.bottom_left.x) r.bottom_left.x = e2->I.bottom_left.x;
    if(e2->I.bottom_left.y < r.bottom_left.y) r.bottom_left.y = e2->I.bottom_left.y;
    if(e2->I.top_right.x > r.top_right.x) r.top_right.x = e2->I.top_right.x;
    if(e2->I.top_right.y > r.top_right.y) r.top_right.y = e2->I.top_right.y;

    return r;
}

// swaps the entries in position pos1 and pos2 in the 'entries' array
void swap_positions(Entry* entries[], int pos1, int pos2){  
    Entry* temp = entries[pos1];
    entries[pos1] = entries[pos2];
    entries[pos2] = temp;
}

/*
From the given M+1 entries in 'entries' array, select two entries which should definitely not be in the same group after splitting.
For every pair, make a rectangle enclosing the two entries of the pair. The pair for which the rectangle has the most free space 
is the desired pair. Store the pointers to the chosen ones in s1 and s2.
*/
void PickSeeds(Entry** s1, Entry** s2, Entry* entries[M+1], int* size){     
    // Assume first two entries are the required ones
    *s1 = entries[0];
    *s2 = entries[1];
    Rectangle r = Enclose_two(*s1, *s2);
    double d = Area(r) - Area((*s1)->I) - Area((*s2)->I);   // Free space in the rectangle enclosing the first two entries

    //Iterate through every possible pair to see if we get a better pair
    for(int i=0; i<M+1; i++){
        for(int j=i+1; j<M+1; j++){
            Rectangle temp_rec = Enclose_two(entries[i], entries[j]);
            double temp_d = Area(temp_rec) - Area(entries[i]->I) - Area(entries[j]->I);
            if(temp_d > d){     // If you get a better pair, update s1 and s2
                *s1 = entries[i];
                *s2 = entries[j];
                d = temp_d;
            }
        }
    }

    // After we got our pair, we remove these two from the 'entries' array by swapping and placing them in the end
    // and decrementing the size variable by 2
    for(int i=0; i<M-1; i++){
        if(entries[i]==*s1){
            swap_positions(entries, i, M-1);
        }
    }
    for(int i=0; i<M-1; i++){
        if(entries[i]==*s2){
            swap_positions(entries, i, M);
        }
    }
    (*size) -= 2;
}

// Expands the given rectangle r just enough to enclose entry e and returns the rectangle
Rectangle Enclose_one(Rectangle r, Entry* e){
    if(e->I.bottom_left.x < r.bottom_left.x) r.bottom_left.x = e->I.bottom_left.x;
    if(e->I.bottom_left.y < r.bottom_left.y) r.bottom_left.y = e->I.bottom_left.y;
    if(e->I.top_right.x > r.top_right.x) r.top_right.x = e->I.top_right.x;
    if(e->I.top_right.y > r.top_right.y) r.top_right.y = e->I.top_right.y;
    return r; 
}

/* For a given array of entries and two groups with MBR r1 and r2, select one entry from the array that has the maximum preference 
for one group over the other    */
Entry* PickNext(Rectangle r1, Rectangle r2, Entry* entries[], int* size){
    if(*size < 1) return NULL;  

    double d1, d2, diff, max_diff;
    Rectangle big_r1, big_r2;
    Entry* temp = entries[0];   // Assume first entry to be the one
    big_r1 = Enclose_one(r1, temp);
    big_r2 = Enclose_one(r2, temp);
    d1 = Area(big_r1) - Area(r1);   // Area expansion required if it gets added to group 1
    d2 = Area(big_r2) - Area(r2);   // Area expansion required if it gets added to group 2

     // The difference between the two expansions will indicate this entry's relative preference for one group over other
    max_diff = (d1-d2 >=0) ? (d1-d2) : (d2-d1); 

    for(int i=1; i < *size; i++){   // Traverse through the rest and see if we get a better one
        big_r1 = Enclose_one(r1, entries[i]);
        big_r2 = Enclose_one(r2, entries[i]);
        d1 = Area(big_r1) - Area(r1);
        d2 = Area(big_r2) - Area(r2);
        diff = (d1-d2 >=0) ? (d1-d2) : (d2-d1);

        if(diff > max_diff){    // If we find a better one, update our temp variable (and max_diff to compare with in other iterations)
            max_diff = diff;
            temp = entries[i];
        }
    }

    // After we get the desired entry, we remove it from the array by placing it at the end and decrementing size variable
    for(int i=0; i< *size-1; i++){
        if(entries[i]==temp){
            swap_positions(entries, i, *size-1);
        }
    }
    --(*size);
    return temp;
}

/*
Given a fully occupied node n1 (containing M entries) and an empty node n2 and another (M+1)th entry e, we need to distribute these 
total of M+1 entries into the two nodes efficiently (inserting each entry into the node whose MBR requires lesser area expansion), 
while also ensuring that both nodes have the minimum required number of entries which is 'm'
*/
void SplitNode(Node* n1, Node* n2, Entry* e){
    Entry* entries[M+1];
    int current_size = M+1;

    for(int i=0; i<n1->occupied; i++) // Store all the M+1 entries in an array 'entries'
        entries[i] = n1->E[i];
    entries[M] = e;

    // Vacate both nodes
    n1->occupied=0;
    n2->occupied=0;
    Rectangle r1, r2, enlarged_r1, enlarged_r2;

    // PickSeeds gives the 1st entry for each of the 2 nodes and stores them in seed1 and seed2 (Also removes these 2 from the 'entries' array)
    Entry *seed1, *seed2;
    PickSeeds(&seed1, &seed2, entries, &current_size);

    add_entry_to_node(seed1, n1);   // adding seed1 to node n1
    r1 = seed1->I;                  // r1 is the MBR for n1
    add_entry_to_node(seed2, n2);   // adding seed2 to node n2
    r2 = seed2->I;                  // r2 is the MBR for n2

    // Keep adding an entry one at a time to one of the nodes until entries array is empty
    while(1){   
        if(current_size<1) return;  // if array is empty, we're done.

        // if n1 has so few entries that all the rest will have to be added to it for it to satisfy the minimum number of entries 'm'
        // then add all of them to n1 and stop (return).
        if(n1->occupied + current_size <= m){
            for(int i=0; i<current_size; i++){
                add_entry_to_node(entries[i], n1);
                return;
            }
        }
        // similarly if n2 has so few entries that all the rest will have to be added to it for it to satisfy the minimum number of entries 
        // 'm' then add all of them to n2 and stop (return).
        else if(n2->occupied + current_size <= m){
            for(int i=0; i<current_size; i++){
                add_entry_to_node(entries[i], n2);
                return;
            }
        }

        //Otherwise pick the next entry from the array that has the most preference for one of the two nodes.
        Entry* temp = PickNext(r1, r2, entries, &current_size);
        enlarged_r1 = r1;
        enlarged_r2 = r2;

        if(temp->I.bottom_left.x < enlarged_r1.bottom_left.x) enlarged_r1.bottom_left.x = temp->I.bottom_left.x;
        if(temp->I.bottom_left.y < enlarged_r1.bottom_left.y) enlarged_r1.bottom_left.y = temp->I.bottom_left.y;
        if(temp->I.top_right.x > enlarged_r1.top_right.x) enlarged_r1.top_right.x = temp->I.top_right.x;
        if(temp->I.top_right.y > enlarged_r1.top_right.y) enlarged_r1.top_right.y = temp->I.top_right.y;

        if(temp->I.bottom_left.x < enlarged_r2.bottom_left.x) enlarged_r2.bottom_left.x = temp->I.bottom_left.x;
        if(temp->I.bottom_left.y < enlarged_r2.bottom_left.y) enlarged_r2.bottom_left.y = temp->I.bottom_left.y;
        if(temp->I.top_right.x > enlarged_r2.top_right.x) enlarged_r2.top_right.x = temp->I.top_right.x;
        if(temp->I.top_right.y > enlarged_r2.top_right.y) enlarged_r2.top_right.y = temp->I.top_right.y;

        // If the entry returned by PickNext has more preference for n1, add it to n1
        if(Area(enlarged_r1)-Area(r1) < Area(enlarged_r2)-Area(r2)){
            add_entry_to_node(temp, n1);
            r1 = enlarged_r1;
        }
        // If it has more preference for n2, add it to n2
        else if(Area(enlarged_r1)-Area(r1) > Area(enlarged_r2)-Area(r2)){
            add_entry_to_node(temp, n2);
            r2 = enlarged_r2;
        }
        else{
            // If both nodes need same area expansion, then add it to the node with smaller current area
            if(Area(r1) < Area(r2)){
                add_entry_to_node(temp, n1);
                r1 = enlarged_r1;
            }
            else if(Area(r1) > Area(r2)){
                add_entry_to_node(temp, n2);
                r2 = enlarged_r2;
            }
            else{
                // If both nodes also have the same current area, add it to the node with fewer entries
                if(n1->occupied < n2->occupied){
                    add_entry_to_node(temp, n1);
                    r1 = enlarged_r1;
                }
                else{
                    add_entry_to_node(temp, n2);
                    r2 = enlarged_r2;
                }
            }
        }
    }
}

// returns the parent node of the node N from the tree rooted at 'root';
// stores the index of the entry in the parent node that points to N in entry_index
Node* Parent(Node* N, Node* root, int* entry_index){
    Node* temp = root;
    if(root->is_leaf) return NULL;  // if the given root itself is a leaf node it cant be the parent of N, hence not found

    for(int i=0; i<temp->occupied; i++){
        if(temp->E[i]->data.p == N){    // Found the parent
            *entry_index = i;
            return temp;
        }
        if(temp->E[i]->data.p){    // Otherwise recursively search in all children of root 
            Node* X = Parent(N, temp->E[i]->data.p, entry_index);
            if(X) return X;
        }
    }
    return NULL;    // if not found
}

// Expands the MBR of entry e so as to accomodate the entries in node n
void Enclose(Entry* e, Node* n){
    int first=1;
    for(int i=0; i<n->occupied; i++){
        Entry* temp_entry = n->E[i];
        if(first){  // For the first entry, just make the MBR of e equal to the MBR of this first entry
            first=0;
            e->I = temp_entry->I;
            continue;
        }
        // Subsequently, expand the MBR of e just enough to enclose the next entry of n
        if(temp_entry->I.bottom_left.x < e->I.bottom_left.x) e->I.bottom_left.x = temp_entry->I.bottom_left.x;
        if(temp_entry->I.bottom_left.y < e->I.bottom_left.y) e->I.bottom_left.y = temp_entry->I.bottom_left.y;
        if(temp_entry->I.top_right.x > e->I.top_right.x) e->I.top_right.x = temp_entry->I.top_right.x;
        if(temp_entry->I.top_right.y > e->I.top_right.y) e->I.top_right.y = temp_entry->I.top_right.y;
    }
}

/*
Whenever an insert is made into a node, AdjustTree is called in order to update the MBRs of the ancestors to accomodate the new entry,
or to propagate node split upwards if needed.
We need to go upward from N and ensure its parent encloses all entries in N (as well as NN if node split led to the formation of NN).
If node split propagation leads to splitting of the root, AdjustTree will return the new sibling of root
*/
Node* AdjustTree(Node* root, Node* N, Node* NN){
    if(N == root){  // if N is root, we have reached the top; just return NN if it is formed due to root split
        return NN;
    }
    int entry_index;
    Node* P = Parent(N, root, &entry_index);
    Entry* EnP = P->E[entry_index];     // EnP is the entry in the parent node P that points to N

    Enclose(EnP, N);    // Update the MBR of EnP to accomodate all entries of N

    Node* PP = NULL;
    if(NN){         // If node splitting led to formation of NN
        Entry* ENN = createEntry(0);    // Create a parent entry ENN to point to NN
        ENN->data.p = NN;
        Enclose(ENN, NN);
        if(P->occupied < M){        // if there is space in the parent node P of N, then add ENN there
            add_entry_to_node(ENN, P);
        }
        else{                       // split parent node P as well to create PP and accomodate ENN
            PP = createNode(0);
            SplitNode(P, PP, ENN);
        }
    }
    return AdjustTree(root, P, PP);    // recurse on the parent(s) to propagate the same thing till the root
}

// Updates the dummy entry Meta_Entry of the tree by expanding its MBR to accomodate all entries of the root node
void resize_mbr_of_tree(R_Tree* the_tree){  
    Enclose(the_tree->Meta_Entry, the_tree->T);
}

// Inserts entry e into theTree
void Insert(Entry* e, R_Tree* theTree){
    Node* L = ChooseLeaf(theTree->T, e);    // chooses the appropriate leaf node to insert into

    Node* LL = NULL;
    if(L->occupied < M){    // if there's space in that node, then add it
        add_entry_to_node(e, L);
    }
    else{                   // otherwise split the node and to produce LL and then accomodate e
        LL = createNode(1);
        SplitNode(L, LL, e);
    }

    // After inserting into a node, call AdjustTree so as to update the ancestors' MBR's to enclose the newly added entry into the node
    Node* second_root = AdjustTree(theTree->T, L, LL);
    if(second_root){    
        // if AdjustTree returned a node, it means node splitting was propagated to the root and the returned node is the new sibling 
        // of root. We make a new root and make the old root and its sibling, its children
        Node* first_root = theTree->T;
        Node* new_root = createNode(0);

        Entry* e1 = createEntry(0);
        e1->data.p = first_root;
        Enclose(e1, first_root);

        Entry* e2 = createEntry(0);
        e2->data.p = second_root;
        Enclose(e2, second_root);

        add_entry_to_node(e1, new_root);
        add_entry_to_node(e2, new_root);
        theTree->T = new_root;
    }
    resize_mbr_of_tree(theTree); // updates the Meta_Entry's rectangle enclosing all entries of root node
}

// Creates a terminal object which is a point in the form of a rectangle with coincident bottom-left and top-right vertices.
Entry* createLeafEntry(double x, double y){
    Entry* temp = createEntry(1);
    Rectangle r;
    r.bottom_left.x = x;
    r.bottom_left.y = y;
    r.top_right.x = x;
    r.top_right.y = y;

    temp->data.obj = r;
    temp->I = r;
    return temp;
}

// Prints a node
void visit(Node* the_node, Rectangle MBR, int level){
    if(the_node->is_leaf){  // if its a leaf, prints all the terminal objects (the 2-D points) stored in it
        fprintf(outp, "\n");
        for(int l=1; l<=level; l++) fprintf(outp, "  ");    // Output indentation
        fprintf(outp, "Printing leaf node: ");
        fprintf(outp, "\n");
        for(int l=1; l<=level; l++) fprintf(outp, "  ");    // Output indentation
        fprintf(outp, "MBR: ");
        fprintf(outp, "\n");
        for(int l=1; l<=level; l++) fprintf(outp, "  ");    // Output indentation
        fprintf(outp, "Bottom-Left: (%ld, %ld)", (long)MBR.bottom_left.x, (long)MBR.bottom_left.y);
        fprintf(outp, "\n");
        for(int l=1; l<=level; l++) fprintf(outp, "  ");    // Output indentation
        fprintf(outp, "Top-Right: (%ld, %ld)", (long)MBR.top_right.x, (long)MBR.top_right.y);
        for(int i=0; i<the_node->occupied; i++){
            fprintf(outp, "\n");
            for(int l=1; l<=level; l++) fprintf(outp, "  ");    // Output indentaion
            fprintf(outp, "2D object %d -- (%ld, %ld)", i+1, (long)the_node->E[i]->data.obj.bottom_left.x, (long)the_node->E[i]->data.obj.bottom_left.y);
        }
    }
    else{   // if its an internal node, prints the MBR of the internal node
        fprintf(outp, "\n");
        for(int l=1; l<=level; l++) fprintf(outp, "  ");    // Output indentation
        fprintf(outp, "Printing internal node: ");
        fprintf(outp, "\n");
        for(int l=1; l<=level; l++) fprintf(outp, "  ");    // Output indentation
        fprintf(outp, "MBR: ");
        fprintf(outp, "\n");
        for(int l=1; l<=level; l++) fprintf(outp, "  ");    // Output indentation
        fprintf(outp, "Bottom-Left: (%ld, %ld)", (long)MBR.bottom_left.x, (long)MBR.bottom_left.y);
        fprintf(outp, "\n");
        for(int l=1; l<=level; l++) fprintf(outp, "  ");    // Output indentation
        fprintf(outp, "Top-Right: (%ld, %ld)", (long)MBR.top_right.x, (long)MBR.top_right.y);
    }
    fprintf(outp, "\n");
}

// Traverses and prints the content of the R tree in preorder
void preorderTraverse(Node* root, Rectangle MBR, int level){
    visit(root, MBR, level);
    if(!(root->is_leaf)){
        for(int i=0; i<root->occupied; i++){
            preorderTraverse(root->E[i]->data.p, root->E[i]->I, level+1);
        }
            
    }
}

int main(){

    R_Tree* myTree = create_R_Tree();

    FILE* fp = fopen("./data.txt", "r");

    double x, y;
    for(int i=1; i<=21 && fscanf(fp, "%lf %lf\n", &x, &y)!=EOF; i++){
        Entry* e = createLeafEntry(x, y);
        Insert(e, myTree);
    }
    fclose(fp);

    outp = fopen("./treecontent.txt", "w");
    preorderTraverse(myTree->T, myTree->Meta_Entry->I, 0);
    fclose(outp);
}