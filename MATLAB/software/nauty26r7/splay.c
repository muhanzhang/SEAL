/* splay.c  - code for splay trees    Version of August 18, 2001.
 * Author: Brendan McKay  bdm@cs.anu.edu.au

   This file is not meant to be compiled separately, but to be
   #included into other programs.  Use it like this:

   1. Define a node type SPLAYNODE.  It must be a structure that
   contains at least the pointer fields left, right and parent of
   type SPLAYNODE*.
   Also define a macro SPLAYNODESIZE giving the size of an object
   of type SPLAYNODE, unless  sizeof(SPLAYNODE) is adequate.

   2. Declare a variable of type SPLAYNODE* to point to the root
   of the tree, and initialise it to NULL.

   3. Declare SCAN_ARGS to be the additional arguments needed for
   splay_scan(), including a leading comma.

   4. Declare ACTION(p) for what splay_scan() should do for node p.

   5. Declare INSERT_ARGS to be the additional arguments needed 
   for splay_insert(), including a leading comma.

   6. Declare COMPARE(p) to compare INSERT_ARGS or LOOKUP_ARGS to the
   contents of node p.  <0, 0, >0 if INSERT_ARGS is greater, equal,
   less, than p.  This has to be an expression with a value, so you will
   need to make it a procedure call if the comparison is complicated.

   If you are using something like strcmp, the correct order is
   strcmp( INSERT_ARGS, p ).

   7. Declare PRESENT(p) for what to do if INSERT_ARGS is already
   present, in node p.  There is a spare int variable i available.
   Typically, this might update some data in the node p.

   8. Declare NOT_PRESENT(p) for what to do if INSERT_ARGS is not
   in the tree and p is a fresh node to hold it.  No need to set
   the left, right, and parent fields.  Use i here too if you like.
   Typically, this might initialise the data in node p.

	PRESENT(p) and NOT_PRESENT(p) should not manipulate the
        tree pointers.  However, each of them can include a
        return if you don't want to change the tree.  In the
        case of NOT_PRESENT(p), do free(p) before returning.

        In the case of PRESENT(p), it is also legal to delete the
        node from the tree using SCAN_DELETE(to_root,p).  In that
        case you MUST return immediately afterwards.

   9. Declare LOOKUP_ARGS to be the additional arguments needed
   for splay_lookup(), including a leading comma.  The default
   for LOOKUP_ARGS is to be the same as INSERT_ARGS.

  10. #include "splay.c"


  Calls:

    Suppose "root" is name of the variable described in step 2.

    There is no need to initialise the tree.  Step 2 did that already.

    To insert something in the tree:
        splay_insert(&root, ...stuff...)
    where "stuff" is the stuff you want to insert, declared as INSERT_ARGS.
    If the key (some part of the stuff decided by you) is present in an
    existing tree node p, PRESENT(p) is executed.  Otherwise, a new tree
    node p is created and NOT_PRESENT(p) is executed.

    To look up something in the tree:
        splay_lookup(&root, ...stuff...)
    where "stuff" is the stuff you want to find, declared as LOOKUP_ARGS.
    It will return a pointer to the tree node with the right key, or NULL
    if there is no such tree node.

    To do something for each node of the tree:
        splay_scan(root, ...stuff...)
    where "stuff" is anything you like (including nothing).  This will
    execute ACTION(p) for each node p in inorder.

    To delete the node p (which MUST be in the tree:
        splay_delete(&root, p)
    Nothing happens if p is NULL, so you can use
        splay_delete(&root, splay_lookup(&root, ...stuff...))
    to delete a node, if any, containing stuff.

  It is possible to have splay trees of several types in the same
  program.  Just include "splay.c" several times, with the procedure
  names SPLAY, SPLAY_SCAN, SPLAY_LOOKUP, SPLAY_INSERT, SPLAY_DELETE
  defined to distinct names.  You have to redefine them all even if
  you aren't using them all.
*/

#define S_A 0
#define S_L 1
#define S_R 2
#define S_LL 3
#define S_LR 4
#define S_RL 5
#define S_RR 6

#ifndef SPLAYNODESIZE
#define SPLAYNODESIZE sizeof(SPLAYNODE)
#endif

#ifndef LOOKUP_ARGS
#define LOOKUP_ARGS INSERT_ARGS
#endif

#ifndef SPLAY
#define SPLAY splay
#define SPLAY_SCAN splay_scan
#define SPLAY_LOOKUP splay_lookup
#define SPLAY_INSERT splay_insert
#define SPLAY_DELETE splay_delete
#endif

/*********************************************************************/

void
SPLAY_SCAN(SPLAYNODE *root SCAN_ARGS)
/* Do ACTION(p) for each node of the tree, in inorder.  Nonrecursive! */
{
    int code;
    SPLAYNODE *p;

    p = root;
    code = S_A;

    while (p)
    {
	switch (code)    /* deliberate flow-ons */
	{
	 case S_A:
	    if (p->left)
	    {
		p = p->left;
		break;
	    }
	 case S_L:
	    ACTION(p);
	    if (p->right)
	    {
		p = p->right;
		code = S_A;
		break;
	    }
	 case S_R:
	    if (p->parent && p->parent->left == p) code = S_L;
	    else                                   code = S_R;
	    p = p->parent;
	    break;
	}
    }
}

/*********************************************************************/

static void
SPLAY(SPLAYNODE *p)
/* Splay the node p.  It becomes the new root. */
{
    SPLAYNODE *q,*r,*s;
    SPLAYNODE *a,*b,*c;
    int code;

#define LCHILD(x,y) {(x)->left = y; if (y) (y)->parent = x;}
#define RCHILD(x,y) {(x)->right = y; if (y) (y)->parent = x;}

    while (p->parent)
    {
	a = p->left;
	b = p->right;
	q = p->parent;
	if (q->left == p)
	{
	    code = S_L;
	    c = q->right;
	}
	else
	{
	    code = S_R;
	    c = q->left;
	}
	r = q->parent;
	if (r)
	{
	    if (r->left == q) code = (code == S_L ? S_LL : S_LR);
	    else              code = (code == S_L ? S_RL : S_RR);
	    s = r->parent;
	    p->parent = s;
	    if (s)
	    {
		if (s->left == r) s->left = p;
		else              s->right = p;
	    }
	}
	else
	{
	    p->parent = NULL;
	}
	
	switch (code)
	{
	 case S_L:
	    RCHILD(p,q);
	    LCHILD(q,b);
	    break;
	 case S_R:
	    LCHILD(p,q);
	    RCHILD(q,a);
	    break;
	 case S_LL:
	    RCHILD(p,q);
	    RCHILD(q,r);
	    LCHILD(q,b);
	    LCHILD(r,c);
	    break;
	 case S_RR:
	    LCHILD(p,q);
	    LCHILD(q,r);
	    RCHILD(r,c);
	    RCHILD(q,a);
	    break;
	 case S_LR:
	    LCHILD(p,q);
	    RCHILD(p,r);
	    RCHILD(q,a);
	    LCHILD(r,b);
	    break;
	 case S_RL:
	    LCHILD(p,r);
	    RCHILD(p,q);
	    RCHILD(r,a);
	    LCHILD(q,b);
	    break;
	}
    }
}

/*********************************************************************/

void
SPLAY_INSERT(SPLAYNODE **to_root  INSERT_ARGS)
/* Do insertion operation.  On return, the object being inserted
   is at the root of the tree regardless of whether a new node
   needed to be created for it. */
{
    int i,cmp;
    SPLAYNODE *p,*ppar,*new_node;

    p = *to_root;
    cmp = 0;

    while (p != NULL)
    {
        cmp = COMPARE(p);
        if (cmp == 0)
        {
	    PRESENT(p);
	    SPLAY(p);
	    *to_root = p;
            return;
        }
        else if (cmp < 0)
        {
            ppar = p;
            p = p->left;
        }
        else
        {
            ppar = p;
            p = p->right;
        }
    }

    if ((new_node = (SPLAYNODE*)malloc(SPLAYNODESIZE)) == NULL)
    {
        fprintf(stderr,">E malloc failed in splay_insert()\n");
        exit(1);
    }

    NOT_PRESENT(new_node);

    new_node->left = new_node->right = NULL;

    if (cmp == 0)
    {
        *to_root = new_node;
	new_node->parent = NULL;
    }
    else if (cmp < 0)
    {
        ppar->left = new_node;
	new_node->parent = ppar;
    }
    else
    {
        ppar->right = new_node;
	new_node->parent = ppar;
    }

    SPLAY(new_node);
    *to_root = new_node;
}

/*********************************************************************/

SPLAYNODE*
SPLAY_LOOKUP(SPLAYNODE **to_root  LOOKUP_ARGS)
/* Do a look-up operation.  If found, return a pointer to the
   node containing it.  If not, return NULL. */
{
    int i,cmp;   /* i is available for COMPARE */
    SPLAYNODE *p;

    p = *to_root;
    cmp = 0;

    while (p != NULL)
    {
        cmp = COMPARE(p);
        if (cmp == 0)
        {
	    SPLAY(p);
	    *to_root = p;
            return p;
        }
        else if (cmp < 0)
            p = p->left;
        else
            p = p->right;
    }

    return NULL;
}

/*********************************************************************/

void
SPLAY_DELETE(SPLAYNODE **to_root, SPLAYNODE *p)
/* Remove node p from the tree and free it. */
{
    SPLAYNODE *q;

    if (p == NULL) return;

    SPLAY(p);
    *to_root = p;

    /* Now we have to delete the root. */

    /* No right child (includes no children). */

    if (!p->right)
    {
	*to_root = p->left;
        if (p->left) p->left->parent = NULL;
	free(p);
	return;
    }

    /* right child but no left child */

    if (!p->left)
    {
        *to_root = p->right;
        p->right->parent = NULL;
        free(p);
        return;
    }  

    /* both children exist */

     for (q = p->left; q->right; q = q->right) {}

     if (q->left) q->left->parent = q->parent;
     if (q->parent == p) q->parent->left = q->left;
     else                q->parent->right = q->left;

     q->left = p->left;
     q->right = p->right;
     q->parent = NULL;
     if (p->left) p->left->parent = q;
     if (p->right) p->right->parent = q;
     *to_root = q;
     free(p);
}

/*********************************************************************/

/* The following shows the tree structure for debugging purposes.
   If you define SPLAY_DUMP you must also define DUMP_ARGS,
   DUMP_LEFT, DUMP_RIGHT and DUMP_ACTION(p). */

#ifdef SPLAY_DUMP
void
SPLAY_DUMP(SPLAYNODE *p DUMP_ARGS)
{
    int i;
    
    if (p == NULL) return;

    if (p->right && p->right->parent != p)
	fprintf(stderr,"parent misaligned at %p-%p ************\n",p,p->right);
    if (p->left && p->left->parent != p)
	fprintf(stderr,"parent misaligned at %p-%p ************\n",p,p->left);

    SPLAY_DUMP(p->right  DUMP_RIGHT);
    DUMP_ACTION(p);
    SPLAY_DUMP(p->left  DUMP_LEFT);
}
#endif
