
# coding: utf-8

# In[114]:


import math


# In[115]:


# First Exercise 
def get_cone_volume(radius, height):
        """The function returns the volume of a cone given, in order, radius and height as input"""
        return (height/3)*math.pi*(radius**2)


# In[69]:


# Calculate the factorial of a number n in a recursive way
def recursive_factorial(n):
    """The function calculates the factorial of a number n in a recursive way. It takes the number n as an input. """
    if n == 0:
        return 1
    else:
        return (n)* recursive_factorial(n-1)
        


# In[70]:


# Calculate the factorial of a number n in an iterative way
def factorial(n):
    """ The function calculates the factorial of a number n in an iterative way. It takes the number n as an input. """
    if n==0:
        return 1
    else: 
        calculated_factorial = 1 
        while n > 0: 
            calculated_factorial*=n
            n-=1
        return(calculated_factorial)


# In[124]:


def count_down(n, odd=False):
    """ This function counts down the number by printing them, from the input n to 0.  The compulsory input is the number from which to count down. There is an optional input, which is a boolean: if set to True the funciton will print only Odd numbers.  """
    while n >=0 : 
        if( odd == True): 
            if (n%2)!= 0:
                print(n)
        else: 
            print(n)
        n-=1


# In[95]:


def get_final_price(price, discount_percentage=10): 
    """Return the final price after applying the discount percentage"""
    return price-( price* discount_percentage / 100)

