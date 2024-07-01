from math import log10

def KL(p, q):
  answer = p * log10(p / q) + (1 - p) * log10((1 - p) / (1 - q))
  return(answer)

def VoI_log_updatable(pu, puc, pc, punotc):
  
  #P(U|c)P(c) and P(U)P(c) terms
  puc_pc = puc * pc * log10((puc * pc )/( pu * pc)) + (1 - puc) * (pc) * log10(((1 - puc) * (pc))/ ((1 - pu) * (pc)))
  
  #P(U|not c)P(not c) and P(U)P(not c) terms
  punotc_pnotc = punotc * (1-pc) * log10((punotc * (1-pc))/ (pu * (1-pc))) + (1 - punotc) * (1-pc) * log10(((1 - punotc) * (1-pc))/ ((1 - pu) * (1-pc)))
  
  answer = puc_pc  + punotc_pnotc
  return answer

def punotc(pc, puc, pu):
  return (pu - puc * pc) / (1 - pc)

def VoI_log(pu, puc, pc):
  punotc = (pu - puc * pc) / (1 - pc)
  l_puc_pu = KL(puc, pu)
  # KL divergence between P(U|¬c) and P(U)
  l_punotc_pu = KL(punotc, pu)
  # EV(KL divergence)
  answer = l_puc_pu * pc + l_punotc_pu * (1 - pc)
  return(answer)