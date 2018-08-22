struct subset{
    int parent;
    int rank;
};

Subset union(Subset s1, Subset s2) {
  if (s1.rank < s2.rank) {
    s1.parent = s2;
  } else if (s1.rank > s2.rank) {
    s2.parent = s1;
  } else {
    s1.parent = s2;
    s1.rank ++;
  }
}
