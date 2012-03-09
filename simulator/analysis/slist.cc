#include "slist.h"

Sllist new_sllist()
{
  Sllist s;
  s= (Sllist)malloc(sizeof(struct sllist));
  s->val=0.0;
  s->next=NULL;
  s->row=0;
  s->col=0;
  return s;
}

Sllist insert_sllist(int nrow, int ncol, double dval, Sllist head)
{
  Sllist prev_ptr=NULL,tmp_ptr=NULL,in_ptr=NULL,new_ptr=NULL;
  in_ptr=NULL;
  int flag=0;
  nrow--;
  ncol--;
  if( head->val==0.0)
  {
		head->row=nrow;
		head->col=ncol;
		head->val = dval;
		return head;
  }
  else
  {
    tmp_ptr=head;
    while(in_ptr==NULL)
    {
      if(tmp_ptr)
      {
        if(tmp_ptr->row<nrow)
        {  // nrow>row
          prev_ptr=tmp_ptr;
          tmp_ptr = tmp_ptr->next;
        }
        else if(tmp_ptr->row==nrow)
        {
          if(tmp_ptr->col==ncol)
          {  //same row and col
            tmp_ptr->val+=dval;
            return head;
          }
          else
          {//same row check for col
            if(tmp_ptr->col<ncol)
            { //ncol>col
              prev_ptr=tmp_ptr;
              tmp_ptr=tmp_ptr->next;
            }
            else 
            { //ncol<col
              if(prev_ptr==NULL)
              {
                flag=1;     //insert before head
                break;
              }
              else
                in_ptr=prev_ptr;
            }
          } //end of col check
        }
        else
        {  //nrow<row
          if(prev_ptr==NULL)
          {
            flag=1;
            break;
          }
          else
            in_ptr=prev_ptr;
        }
      }
      else in_ptr=prev_ptr;
    }
    new_ptr = (Sllist)malloc(sizeof(struct sllist));
    if(flag)
    { //insert before head
      new_ptr->next=head;
      head=new_ptr;
    }
    else
    {
      new_ptr->next=in_ptr->next;
      in_ptr->next=new_ptr;
    }
    new_ptr->row=nrow;
    new_ptr->col=ncol;
    new_ptr->val=dval;
    return head;
  }
}

void sllist_print(Sllist head)
{
  Sllist tmp_ptr;
  tmp_ptr=head;
  while(tmp_ptr)
  {
    cout<<"Mp("<<tmp_ptr->row<<","<<tmp_ptr->col<<") "<<tmp_ptr->val<<endl;
    tmp_ptr=tmp_ptr->next;
  }
}
