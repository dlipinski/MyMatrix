# MyMatrix
Implementation for matrix, matrix operations and special data type (fraction)
### Operacje na macierzach
Zdefiniowany został parametryzowany szablon klasy C++ reprezentujący macierz nad
ciałem liczb rzeczywistych oraz stworzone przeciążone operatory dodawania, odejmowania,
mnożenia i dzielenia.
### Operacje na ułamkach
Zdefiniowany został szablon klacy C++ reprezentujący ułamek, wraz z przeciążonymi
operatorami, składający się z licznika I mianownika, które to są zmiennymi typu BigInteger
(zewnętrzna klasa).
### Testy poprawności
Wszystkie testy (nie licząc metody Gaussa) przeprowadzone zostały używając
następujących typów reprezentujących
liczbę rzeczywistą:

• typu pojedynczej precyzji: float,

• typu podwójnej precyzji: double

• typu własnego MyFraction, który przechowuje liczbę w postaci ułamka liczb typu BigInteger.
### Dodawanie, mnożenie oraz implementacja metody Gaussa
Dla losowych macierzy kwadratowych A, B, C I wektora X wykonane zostały testy badające
poprawność (błędy) I wydajność (czas działania następujących operacji):

• A * X,

• (A + B + C) * X,

• A * (B * C),

• metoda Gaussa z pełnym wyborem elementu podstawowego,

• metoda Gaussa z częściowym wyborem elementu podstawowego.

Wyniki i czas działania powyższych operacji zostały porównane z wynikami uzyskanymi
przu użyciu klasy Matrix z biblioteki Eigen3.
