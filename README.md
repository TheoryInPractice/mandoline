# Mandoline

To install, run

```
pip install git+https://github.com/TheoryInPractice/mandoline
```

To count a pattern graph, mandoline needs to first generate a so-called 'counting dag' from the pattern. For example, to generate a pattern graph of the 'cricket' graph, run

```
  mandoline decompose example-graphs/cricket.txt --output cricket.dag
```

Now we can count the number of times this graph appears in the `codeminar` graph by running 

```
  mandoline count cricket.dag example-graphs/codeminer.txt.gz
```
