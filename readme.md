# 2020.04.10

    The code will not compile on the iMac that I am using. I don't wish the
    quagmire to refactor all the code, but I do wish to clean it up
    substantially. For now, I am moving file and classes from names like `io*`
    to `tu*`, with the general idea to take what I need in a clean manner as I 
    need it.



* Two types of code:



    1. `src/io*`
        Contain all code that is compiled by ROOT and can be used interactively in ROOT
        from the ./lib directory
    2. `src/oi*`
        Contain code that is not compiled directly by ROOT, and is used to make object files
        under ./obj directory
