#ifndef __HEAP_H
#define __HEAP_H

#include <vector>




class CHeapItem
{
  public:
    float Value;
    int Index; 

  private:
};




class CHeap
{
  public:
    std::vector<CHeapItem> Items;
    std::vector<int> Indices;
    int CurrentNumItems;

    void Put(const int Index, const float Value);
    void Update(const int Index, const float Value);
		void Clear();
    void Initialize(const int NumItems);

    // not used 
    void Remove(const int Index);
    int Get();

    CHeap()
    {
      this->CurrentNumItems = 0;
    }



  private:
	  inline void SwapItems(const int i, const int j);
    void UpHeap(const int i);
    void DownHeap(const int i);
};




inline void CHeap::SwapItems(const int i, const int j)
{
  CHeapItem Temp;

  Temp = this->Items[i];
  this->Items[i] = this->Items[j];
  this->Items[j] = Temp;

  this->Indices[this->Items[i].Index] = i;
  this->Indices[this->Items[j].Index] = j;
}



#endif
