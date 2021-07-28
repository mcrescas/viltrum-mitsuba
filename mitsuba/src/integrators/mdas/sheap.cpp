#include "heap.h"





void CHeap::Clear()
{
	this->Items.clear();
	this->Indices.clear();
}






void CHeap::UpHeap(const int i)
{
  int k = i;

  while (k != 0)
  {
    int j = (k - 1) / 2;

    if (this->Items[k].Value <= this->Items[j].Value)
    {
      break;
    }

    this->SwapItems(k, j);

    k = j;
  }
}





void CHeap::DownHeap(const int i)
{
  int j = (i * 2 + 1);
  int k = i;
 
	int m = this->CurrentNumItems;
  while (j < m)
  {
    if ((j + 1) < m) 
		{
			if (this->Items[j + 1].Value > this->Items[j].Value)
			{
				j++;
			}
		}

    if (this->Items[j].Value <= this->Items[k].Value)
    {
      break;
    }

    this->SwapItems(k, j);

    k = j;
    j = (k * 2 + 1);
  }
}





void CHeap::Put(const int Index, const float Value)
{
	CHeapItem NewItem;
	NewItem.Index = Index;
	NewItem.Value = Value;
	//this->Items.push_back(NewItem);
	this->Items[this->CurrentNumItems] = NewItem;
	this->Indices[Index] = (this->CurrentNumItems);
  //if (this->Indices.size() <= Index)
  //{
  //  this->Indices.resize(Index + 1);
  //}
  this->UpHeap(this->CurrentNumItems);

  this->CurrentNumItems++;
}





int CHeap::Get()
{
  int Index = this->Items[0].Index;

  this->Items[0] = this->Items[this->Items.size() - 1];
  this->Indices[this->Items[0].Index] = 0;
	this->Items.pop_back();
	//this->Indices.pop_back();

	this->DownHeap(0);

  return Index;
}





void CHeap::Update(const int Index, const float Value)
{
  int i = this->Indices[Index];
  float OldValue = this->Items[i].Value;
  this->Items[i].Value = Value;

  if (Value > OldValue)
  {
    this->UpHeap(i);
  }
  else
  {
    this->DownHeap(i);
  }
}





void CHeap::Remove(const int Index)
{
  int i = this->Indices[Index];
  float OldValue = this->Items[i].Value;
  this->Items[i] = this->Items[this->Items.size() - 1];
  this->Indices[this->Items[i].Index] = i;
	this->Items.pop_back();

  if (Items[i].Value > OldValue)
  {
    this->UpHeap(i);
  }
  else
  {
    this->DownHeap(i);
  }
}
