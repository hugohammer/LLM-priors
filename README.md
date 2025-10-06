# LLM-priors

This repository contains the code associated with the paper:  

[**Using Large Language Models to Suggest Informative Prior Distributions in Bayesian Regression Analysis**](https://www.nature.com/articles/s41598-025-18425-9)  
*Michael A. Riegler, Kristoffer H. Hellton, Vajira Thambawita, and Hugo L. Hammer*  
Published in *Scientific Reports*, 2025.

---

## Summary

This paper explores the potential of using large language models (LLMs) to select informative prior distributions in Bayesian regression analysis.  

We evaluated three different LLMs ChatGPT, Gemini, and Claude and found that the models were able to correctly identify the direction of associations between predictors and the response variable. This demonstrates great potential for leveraging LLMs to select informative priors. However, a significant challenge remains in calibrating the width of these priors. The LLMs showed tendencies toward both overconfidence and underconfidence.

---

## Citation

If you use this code or paper in your research, please cite it as:

```bibtex
@article{riegler2025using,
  title={Using large language models to suggest informative prior distributions in Bayesian regression analysis},
  author={Riegler, Michael A and Hellton, Kristoffer H and Thambawita, Vajira and Hammer, Hugo L},
  journal={Scientific Reports},
  volume={15},
  number={1},
  pages={33386},
  year={2025},
  publisher={Nature Publishing Group UK London}
}
