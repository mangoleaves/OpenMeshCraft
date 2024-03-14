Rem FOR %%f in (.\direct\*.txt) DO ..\..\..\build\bin\Release\PredicatesGenerator.exe %%f -a -ft -it -et -ex
FOR %%f in (.\indirect\*.txt) DO ..\..\..\build\bin\Release\PredicatesGenerator.exe %%f -a -ft -it -et -ex
FOR %%f in (.\lambda\*.txt) DO ..\..\..\build\bin\Release\PredicatesGenerator.exe %%f -a -ft -it -et -ex
COPY indirect_predicates.h .\Generated
COPY indirect_predicates.hpp .\Generated
DEL indirect_predicates.h
DEL indirect_predicates.hpp
