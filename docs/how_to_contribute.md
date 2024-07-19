# How to contribute
## Contributing code
### Coding standards
- Update VERSION value in pattools/INFO.yaml if the code is changed 
- When a milestone is reached, such as a new release, use `git tag` to mark specific commits.
```
git tag -a v0.1.16 -m "Release version 0.1.16"
```

- Try to use the basic code provided in the project to perform common operations, such as cutting pat files, 
   reading pat files, reading tabix index files, and outputting files. This approach is beneficial for the 
   overall optimization of the project.
